#include "PitchScaling.h"
#include <cmath>
#include <complex>

namespace audio {

PitchScaler::PitchScaler(const EngineConfig& config)
    : fftSize(config.fftSize)
    , hopSize(config.hopSize)
    , pitchFactor(std::pow(2.0f, config.pitchShift / 12.0f))
    , sampleRate(config.sampleRate)
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

            // Initialize phase accumulation buffers
            lastPhase.resize(fftSize / 2 + 1);
            sumPhase.resize(fftSize / 2 + 1);
            std::fill(lastPhase.begin(), lastPhase.end(), 0.0f);
            std::fill(sumPhase.begin(), sumPhase.end(), 0.0f);
        }

PitchScaler::~PitchScaler() {
    fftw_destroy_plan(fftPlan);
    fftw_destroy_plan(ifftPlan);
    fftw_free(fftIn);
    fftw_free(fftOut);
    fftw_free(ifftIn);
    fftw_free(ifftOut);
}

void PitchScaler::processFrame(const std::vector<float>& input,
                              std::vector<float>& output,
                              size_t inputOffset,
                              size_t outputOffset)
{
            // Apply window and perform FFT
            for (size_t i = 0; i < fftSize; ++i) {
                fftIn[i] = input[inputOffset + i] * window[i];
            }
            fftw_execute(fftPlan);

            // Convert to magnitude/phase representation
            Eigen::VectorXf magnitude(fftSize / 2 + 1);
            Eigen::VectorXf phase(fftSize / 2 + 1);
            Eigen::VectorXf frequency(fftSize / 2 + 1);

            for (size_t bin = 0; bin <= fftSize/2; ++bin) {
                const float re = fftOut[bin][0];
                const float im = fftOut[bin][1];
                
                magnitude(bin) = std::sqrt(re * re + im * im);
                phase(bin) = std::atan2(im, re);

                // Calculate true frequency using phase difference
                float phaseDiff = phase(bin) - lastPhase[bin];
                // Wrap phase difference to [-π, π]
                while (phaseDiff > M_PI) phaseDiff -= 2.0f * M_PI;
                while (phaseDiff < -M_PI) phaseDiff += 2.0f * M_PI;

                const float expectedPhase = 2.0f * M_PI * bin * hopSize / fftSize;
                const float deviation = phaseDiff - expectedPhase;
                
                frequency(bin) = (expectedPhase + deviation) * sampleRate / (2.0f * M_PI * hopSize);
                lastPhase[bin] = phase(bin);
            }

            // Modify frequency/magnitude for pitch shifting
            for (size_t bin = 0; bin <= fftSize/2; ++bin) {
                const float pitchedFreq = frequency(bin) * pitchFactor;
                // Removed unused variable pitchedBin
                
                // Accumulate phase
                sumPhase[bin] += 2.0f * M_PI * pitchedFreq * hopSize / sampleRate;

                // Calculate new complex values
                const float newPhase = sumPhase[bin];
                ifftIn[bin][0] = magnitude(bin) * std::cos(newPhase);
                ifftIn[bin][1] = magnitude(bin) * std::sin(newPhase);
            }

            // Apply formant preservation
            if (std::abs(pitchFactor - 1.0f) > 0.01f) {
                preserveFormants(magnitude);
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

void PitchScaler::preserveFormants(const Eigen::VectorXf& originalMagnitude) {
            const int filterOrder = 30;  // LPC order for formant estimation
            
            // Extract spectral envelope using LPC
            Eigen::VectorXf r(filterOrder + 1);  // Autocorrelation coefficients
            for (int i = 0; i <= filterOrder; ++i) {
                float sum = 0.0f;
                for (size_t j = 0; j < fftSize/2 - i; ++j) {
                    sum += originalMagnitude(j) * originalMagnitude(j + i);
                }
                r(i) = sum;
            }

            // Levinson-Durbin recursion to get LPC coefficients
            Eigen::VectorXf a(filterOrder + 1);
            a(0) = 1.0f;
            float e = r(0);

            for (int i = 1; i <= filterOrder; ++i) {
                float k = r(i);
                for (int j = 1; j < i; ++j) {
                    k -= a(j) * r(i - j);
                }
                k /= e;
                
                a(i) = k;
                for (int j = 1; j < i; ++j) {
                    // Removed unused variable aj and use a(j) directly
                    a(j) -= k * a(i - j);
                }
                e *= (1.0f - k * k);
            }

            // Apply formant correction
            for (size_t bin = 0; bin <= fftSize/2; ++bin) {
                const float freq = bin * sampleRate / fftSize;
                const float pitchedFreq = freq * pitchFactor;
                
                // Evaluate formant envelope at original and pitched frequencies
                std::complex<float> originalResponse(1.0f, 0.0f);
                std::complex<float> pitchedResponse(1.0f, 0.0f);
                
                for (int i = 1; i <= filterOrder; ++i) {
                    const float originalPhase = -2.0f * M_PI * freq * i / sampleRate;
                    const float pitchedPhase = -2.0f * M_PI * pitchedFreq * i / sampleRate;
                    
                    originalResponse += a(i) * std::polar(1.0f, originalPhase);
                    pitchedResponse += a(i) * std::polar(1.0f, pitchedPhase);
                }

                // Apply correction factor
                const float correctionFactor = std::abs(originalResponse / pitchedResponse);
                ifftIn[bin][0] *= correctionFactor;
                ifftIn[bin][1] *= correctionFactor;
            }
        }


PitchScalingProcessor::PitchScalingProcessor(const EngineConfig& config)
    : scaler(std::make_unique<PitchScaler>(config))
    , hopSize(config.hopSize)
    , fftSize(config.fftSize)
{}

void PitchScalingProcessor::process(const std::vector<float>& input, std::vector<float>& output) {
    const size_t numFrames = 1 + (input.size() - fftSize) / hopSize;
    
    // Initialize output buffer
    output.resize(input.size());
    std::fill(output.begin(), output.end(), 0.0f);

    // Process frame by frame
    for (size_t frame = 0; frame < numFrames; ++frame) {
        const size_t inputOffset = frame * hopSize;
        const size_t outputOffset = inputOffset;  // No time stretching
        
        scaler->processFrame(input, output, inputOffset, outputOffset);
    }
}

} // namespace audio