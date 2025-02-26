#include "../../include/stretch_engine.h"
#include "STNDecomposition.h"
#include <algorithm>
#include <cmath>

namespace audio {
namespace detail {

STNDecompositionImpl::STNDecompositionImpl(const EngineConfig& config)
    : fftSize(config.fftSize)
    , hopSize(config.hopSize)
    , sinesMagnitude(fftSize / 2 + 1)
    , sinesPhase(fftSize / 2 + 1)
    , transientsMagnitude(fftSize / 2 + 1)
    , transientsPhase(fftSize / 2 + 1)
    , noiseMagnitude(fftSize / 2 + 1)
    , noisePhase(fftSize / 2 + 1)
{
    if (fftSize == 0 || hopSize == 0 || hopSize > fftSize) {
        throw std::runtime_error("Invalid FFT or hop size configuration");
    }

    // Initialize Eigen vectors
    sinesMagnitude.setZero();
    sinesPhase.setZero();
    transientsMagnitude.setZero();
    transientsPhase.setZero();
    noiseMagnitude.setZero();
    noisePhase.setZero();

    // Allocate FFTW plans and buffers
    fftIn = fftw_alloc_real(fftSize);
    if (!fftIn) {
        throw std::runtime_error("Failed to allocate FFT input buffer");
    }

    fftOut = fftw_alloc_complex(fftSize / 2 + 1);
    if (!fftOut) {
        fftw_free(fftIn);
        throw std::runtime_error("Failed to allocate FFT output buffer");
    }

    ifftIn = fftw_alloc_complex(fftSize / 2 + 1);
    if (!ifftIn) {
        fftw_free(fftOut);
        fftw_free(fftIn);
        throw std::runtime_error("Failed to allocate IFFT input buffer");
    }

    ifftOut = fftw_alloc_real(fftSize);
    if (!ifftOut) {
        fftw_free(ifftIn);
        fftw_free(fftOut);
        fftw_free(fftIn);
        throw std::runtime_error("Failed to allocate IFFT output buffer");
    }

    fftPlan = fftw_plan_dft_r2c_1d(fftSize, fftIn, fftOut, FFTW_MEASURE);
    if (!fftPlan) {
        fftw_free(ifftOut);
        fftw_free(ifftIn);
        fftw_free(fftOut);
        fftw_free(fftIn);
        throw std::runtime_error("Failed to create FFT plan");
    }

    ifftPlan = fftw_plan_dft_c2r_1d(fftSize, ifftIn, ifftOut, FFTW_MEASURE);
    if (!ifftPlan) {
        fftw_destroy_plan(fftPlan);
        fftw_free(ifftOut);
        fftw_free(ifftIn);
        fftw_free(fftOut);
        fftw_free(fftIn);
        throw std::runtime_error("Failed to create IFFT plan");
    }
    
    // Initialize analysis windows
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
}

STNDecompositionImpl::~STNDecompositionImpl() {
    if (ifftPlan) fftw_destroy_plan(ifftPlan);
    if (fftPlan) fftw_destroy_plan(fftPlan);
    if (ifftOut) fftw_free(ifftOut);
    if (ifftIn) fftw_free(ifftIn);
    if (fftOut) fftw_free(fftOut);
    if (fftIn) fftw_free(fftIn);
}

void STNDecompositionImpl::decompose(
    const std::vector<float>& input,
    std::vector<float>& sines,
    std::vector<float>& transients,
    std::vector<float>& noise)
{
    if (input.size() < fftSize) {
        throw std::runtime_error("Input buffer size must be at least FFT size");
    }

    const size_t numFrames = 1 + (input.size() - fftSize) / hopSize;
    if (numFrames == 0) {
        throw std::runtime_error("Invalid input size or hop size configuration");
    }
    
    // Resize output vectors
    sines.resize(input.size());
    transients.resize(input.size());
    noise.resize(input.size());
    std::fill(sines.begin(), sines.end(), 0.0f);
    std::fill(transients.begin(), transients.end(), 0.0f);
    std::fill(noise.begin(), noise.end(), 0.0f);

    // Temporary vectors for spectral processing
    Eigen::VectorXf magnitude(fftSize / 2 + 1);
    Eigen::VectorXf phase(fftSize / 2 + 1);
    
    // Process frame by frame
    for (size_t frame = 0; frame < numFrames; ++frame) {
        const size_t offset = frame * hopSize;
        
        // Apply window and perform FFT
        for (size_t i = 0; i < fftSize; ++i) {
            fftIn[i] = input[offset + i] * window[i];
        }
        fftw_execute(fftPlan);

        // Convert to magnitude/phase representation
        for (size_t bin = 0; bin <= fftSize/2; ++bin) {
            const float re = fftOut[bin][0];
            const float im = fftOut[bin][1];
            magnitude(bin) = std::sqrt(re * re + im * im);
            phase(bin) = std::atan2(im, re);
        }

        // Perform spectral decomposition
        decomposeSpectrum(magnitude, phase, offset, sines, transients, noise);
    }
}

// Move constructor implementation
STNDecompositionImpl::STNDecompositionImpl(STNDecompositionImpl&& other) noexcept
    : fftSize(other.fftSize)
    , hopSize(other.hopSize)
    , window(std::move(other.window))
    , fftIn(other.fftIn)
    , fftOut(other.fftOut)
    , ifftIn(other.ifftIn)
    , ifftOut(other.ifftOut)
    , fftPlan(other.fftPlan)
    , ifftPlan(other.ifftPlan)
    , sinesMagnitude(std::move(other.sinesMagnitude))
    , sinesPhase(std::move(other.sinesPhase))
    , transientsMagnitude(std::move(other.transientsMagnitude))
    , transientsPhase(std::move(other.transientsPhase))
    , noiseMagnitude(std::move(other.noiseMagnitude))
    , noisePhase(std::move(other.noisePhase))
{
    // Clear other's pointers to prevent double free
    other.fftIn = nullptr;
    other.fftOut = nullptr;
    other.ifftIn = nullptr;
    other.ifftOut = nullptr;
    other.fftPlan = nullptr;
    other.ifftPlan = nullptr;
}

// Move assignment implementation
STNDecompositionImpl& STNDecompositionImpl::operator=(STNDecompositionImpl&& other) noexcept {
    if (this != &other) {
        // Clean up existing resources
        if (ifftPlan) fftw_destroy_plan(ifftPlan);
        if (fftPlan) fftw_destroy_plan(fftPlan);
        if (ifftOut) fftw_free(ifftOut);
        if (ifftIn) fftw_free(ifftIn);
        if (fftOut) fftw_free(fftOut);
        if (fftIn) fftw_free(fftIn);

        // Move resources from other
        window = std::move(other.window);
        fftIn = other.fftIn;
        fftOut = other.fftOut;
        ifftIn = other.ifftIn;
        ifftOut = other.ifftOut;
        fftPlan = other.fftPlan;
        ifftPlan = other.ifftPlan;
        sinesMagnitude = std::move(other.sinesMagnitude);
        sinesPhase = std::move(other.sinesPhase);
        transientsMagnitude = std::move(other.transientsMagnitude);
        transientsPhase = std::move(other.transientsPhase);
        noiseMagnitude = std::move(other.noiseMagnitude);
        noisePhase = std::move(other.noisePhase);

        // Clear other's pointers
        other.fftIn = nullptr;
        other.fftOut = nullptr;
        other.ifftIn = nullptr;
        other.ifftOut = nullptr;
        other.fftPlan = nullptr;
        other.ifftPlan = nullptr;
    }
    return *this;
}

void STNDecompositionImpl::decomposeSpectrum(
    const Eigen::VectorXf& magnitude,
    const Eigen::VectorXf& phase,
    size_t frameOffset,
    std::vector<float>& sines,
    std::vector<float>& transients,
    std::vector<float>& noise)
{
    // Peak detection parameters
    const float peakThreshold = 0.1f;  // Relative threshold for peak detection
    const size_t minPeakDistance = 3;  // Minimum distance between peaks in bins
    
    // Spectral mask buffers
    Eigen::VectorXf sineMask = Eigen::VectorXf::Zero(fftSize/2 + 1);
    Eigen::VectorXf transientMask = Eigen::VectorXf::Zero(fftSize/2 + 1);
    Eigen::VectorXf noiseMask = Eigen::VectorXf::Zero(fftSize/2 + 1);

    // Find spectral peaks for sinusoidal components
    for (size_t bin = 1; bin < fftSize/2; ++bin) {
        if (magnitude(bin) > magnitude(bin-1) && 
            magnitude(bin) > magnitude(bin+1) &&
            magnitude(bin) > peakThreshold * magnitude.maxCoeff()) {
            // Mark as sinusoidal component
            sineMask(bin) = 1.0f;
            
            // Apply peak width
            for (size_t i = 1; i < minPeakDistance && (bin+i) <= fftSize/2; ++i) {
                sineMask(bin+i) = std::max(0.0f, 1.0f - float(i)/minPeakDistance);
            }
            for (size_t i = 1; i < minPeakDistance && (bin-i) > 0; ++i) {
                sineMask(bin-i) = std::max(0.0f, 1.0f - float(i)/minPeakDistance);
            }
        }
    }

    // Detect transients using phase deviation
    Eigen::VectorXf phaseDeviation = Eigen::VectorXf::Zero(fftSize/2 + 1);
    for (size_t bin = 1; bin < fftSize/2; ++bin) {
        float expectedPhase = 2.0f * M_PI * bin * hopSize / fftSize;
        float deviation = std::abs(phase(bin) - expectedPhase);
        while (deviation > M_PI) deviation -= 2.0f * M_PI;
        phaseDeviation(bin) = std::abs(deviation);
        
        if (phaseDeviation(bin) > 0.8f * M_PI) {
            transientMask(bin) = 1.0f;
        }
    }

    // Remaining energy is considered noise
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        float total = sineMask(bin) + transientMask(bin);
        noiseMask(bin) = std::max(0.0f, 1.0f - total);
    }

    // Apply masks to reconstruct components
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        // Convert magnitude/phase back to complex for each component
        const std::complex<double> value(
            magnitude(bin) * std::cos(phase(bin)),
            magnitude(bin) * std::sin(phase(bin))
        );

        // Apply masks to separate components
        if (bin == 0 || bin == fftSize/2) {
            // Handle DC and Nyquist bins (real-only)
            fftOut[bin][0] = value.real();
            fftOut[bin][1] = 0.0;
        } else {
            // Regular bins
            fftOut[bin][0] = value.real();
            fftOut[bin][1] = value.imag();
            // Mirror for negative frequencies (only if not at DC or Nyquist)
            if (bin > 0 && bin < fftSize/2) {
                const size_t mirrorBin = fftSize/2 + 1 - bin;
                if (mirrorBin < fftSize/2 + 1) {
                    fftOut[mirrorBin][0] = value.real();
                    fftOut[mirrorBin][1] = -value.imag();
                }
            }
        }

        // Process sinusoidal component
        ifftIn[bin][0] = fftOut[bin][0] * sineMask(bin);
        ifftIn[bin][1] = fftOut[bin][1] * sineMask(bin);
    }
    
    // Inverse FFT for sinusoidal component
    fftw_execute(ifftPlan);
    for (size_t i = 0; i < fftSize; ++i) {
        const size_t outIndex = frameOffset + i;
        if (outIndex < sines.size()) {
            sines[outIndex] += ifftOut[i] * window[i] / fftSize;
        }
    }

    // Process transient component
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        ifftIn[bin][0] = fftOut[bin][0] * transientMask(bin);
        ifftIn[bin][1] = fftOut[bin][1] * transientMask(bin);
    }
    
    // Inverse FFT for transient component
    fftw_execute(ifftPlan);
    for (size_t i = 0; i < fftSize; ++i) {
        const size_t outIndex = frameOffset + i;
        if (outIndex < transients.size()) {
            transients[outIndex] += ifftOut[i] * window[i] / fftSize;
        }
    }

    // Process noise component
    for (size_t bin = 0; bin <= fftSize/2; ++bin) {
        ifftIn[bin][0] = fftOut[bin][0] * noiseMask(bin);
        ifftIn[bin][1] = fftOut[bin][1] * noiseMask(bin);
    }
    
    // Inverse FFT for noise component
    fftw_execute(ifftPlan);
    for (size_t i = 0; i < fftSize; ++i) {
        const size_t outIndex = frameOffset + i;
        if (outIndex < noise.size()) {
            noise[outIndex] += ifftOut[i] * window[i] / fftSize;
        }
    }
}

} // namespace detail

// Public interface implementation using the Impl
STNDecomposition::STNDecomposition(const EngineConfig& config)
    : pImpl(new detail::STNDecompositionImpl(config)) {}

STNDecomposition::~STNDecomposition() = default;

void STNDecomposition::decompose(const std::vector<float>& input,
                               std::vector<float>& sines,
                               std::vector<float>& transients,
                               std::vector<float>& noise)
{
    pImpl->decompose(input, sines, transients, noise);
}

} // namespace audio
