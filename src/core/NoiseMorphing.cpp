#include "NoiseMorphing.h"
#include <random>
#include <cmath>

namespace audio {

NoiseMorpher::NoiseMorpher(const EngineConfig& config)
    : fftSize(config.fftSize)
    , hopSize(config.hopSize)
    , stretchFactor(config.stretchFactor)
{
    // Initialize FFTW
    fftIn = fftw_alloc_real(fftSize);
    fftOut = fftw_alloc_complex(fftSize / 2 + 1);
    ifftIn = fftw_alloc_complex(fftSize / 2 + 1);
    ifftOut = fftw_alloc_real(fftSize);
    
    fftPlan = fftw_plan_dft_r2c_1d(fftSize, fftIn, fftOut, FFTW_MEASURE);
    ifftPlan = fftw_plan_dft_c2r_1d(fftSize, ifftIn, ifftOut, FFTW_MEASURE);

    // Initialize window
    window.resize(fftSize);
    float sum = 0.0f;
    for (size_t i = 0; i < fftSize; ++i) {
        window[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (fftSize - 1)));
        sum += window[i];
    }
    // Normalize window for perfect reconstruction
    float normFactor = 2.0f * hopSize / sum;
    for (float& w : window) {
        w *= normFactor;
    }

    // Initialize random number generator
    rng.seed(std::random_device{}());
}

NoiseMorpher::~NoiseMorpher() {
    fftw_destroy_plan(fftPlan);
    fftw_destroy_plan(ifftPlan);
    fftw_free(fftIn);
    fftw_free(fftOut);
    fftw_free(ifftIn);
    fftw_free(ifftOut);
}

void NoiseMorpher::processFrame(const std::vector<float>& input,
                               std::vector<float>& output,
                               size_t inputOffset,
                               size_t outputOffset)
{
    // Apply window and perform FFT
    for (size_t i = 0; i < fftSize; ++i) {
        fftIn[i] = input[inputOffset + i] * window[i];
    }
    fftw_execute(fftPlan);

    // Convert to log-magnitude spectrum
    Eigen::VectorXf logMagnitude(fftSize / 2 + 1);
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        const float re = fftOut[bin][0];
        const float im = fftOut[bin][1];
        const float mag = std::sqrt(re * re + im * im);
        logMagnitude(bin) = mag > 1e-10f ? std::log(mag) : -23.0f; // -23 ~= log(1e-10)
    }

    // Interpolate log-magnitude spectra for time-stretching
    Eigen::VectorXf stretchedLogMag = interpolateSpectrum(logMagnitude);

    // Generate noise with stretched spectral envelope
    synthesizeNoise(stretchedLogMag, outputOffset, output);
}

Eigen::VectorXf NoiseMorpher::interpolateSpectrum(const Eigen::VectorXf& logMag) const {
    // For now, implement simple linear interpolation
    // TODO: Implement more sophisticated interpolation methods
    return logMag;
}

void NoiseMorpher::synthesizeNoise(const Eigen::VectorXf& logMag,
                                  size_t outputOffset,
                                  std::vector<float>& output) {
    // Generate white noise with desired spectral envelope
    std::uniform_real_distribution<float> phaseDist(-M_PI, M_PI);
    
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        const float magnitude = std::exp(logMag(bin));
        const float phase = phaseDist(rng);
        
        ifftIn[bin][0] = magnitude * std::cos(phase);
        ifftIn[bin][1] = magnitude * std::sin(phase);
    }

    // Inverse FFT
    fftw_execute(ifftPlan);

    // Apply window and overlap-add
    const float scale = 1.0f / fftSize;
    for (size_t i = 0; i < fftSize; ++i) {
        const size_t outIndex = outputOffset + i;
        const float value = ifftOut[i] * window[i] * scale;

        if (outIndex < output.size()) {
            output[outIndex] += value;
        }
    }
}

NoiseMorphingProcessor::NoiseMorphingProcessor(const EngineConfig& cfg)
    : morpher(std::make_unique<NoiseMorpher>(cfg))
    , hopSize(cfg.hopSize)
    , fftSize(cfg.fftSize)
    , config(cfg)
{}

void NoiseMorphingProcessor::process(const std::vector<float>& input, std::vector<float>& output) {
    // Calculate number of frames
    const size_t numFrames = 1 + (input.size() - fftSize) / hopSize;
    
    // Resize output buffer with overlap
    const size_t outputLength = size_t(input.size() * config.stretchFactor);
    output.resize(outputLength);
    std::fill(output.begin(), output.end(), 0.0f);

    // Process frame by frame
    for (size_t frame = 0; frame < numFrames; ++frame) {
        const size_t inputOffset = frame * hopSize;
        const size_t outputOffset = size_t(inputOffset * config.stretchFactor);
        
        morpher->processFrame(input, output, inputOffset, outputOffset);
    }
}

} // namespace audio
