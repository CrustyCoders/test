#include "../../include/stretch_engine.h"
#include "STNDecomposition.h"
#include "PitchScaling.h"
#include "NoiseMorphing.h"
#include <algorithm>
#include <memory>
#include <mutex>

namespace audio {

class StretchEngine::Impl {
public:
    explicit Impl(const EngineConfig& config)
        : config(config)
        , stnDecomp(config)
        , latency(config.fftSize)
    {
        // Initialize intermediate buffers with proper size
        const size_t bufferSize = config.blockSize * 2;
        sineComponent.resize(bufferSize, 0.0f);
        transientComponent.resize(bufferSize, 0.0f);
        noiseComponent.resize(bufferSize, 0.0f);
        processedSines.resize(bufferSize, 0.0f);
        processedTransients.resize(bufferSize, 0.0f);
        processedNoise.resize(bufferSize, 0.0f);
    }

    void updateConfig(const EngineConfig& newConfig) {
        std::lock_guard<std::mutex> lock(configMutex);
        config = newConfig;
        
        // Create new instance and swap
        STNDecomposition newDecomp(newConfig);
        stnDecomp = std::move(newDecomp);
        configChanged = true;  // Flag to recreate other components
        
        // Update latency
        latency = newConfig.fftSize;
    }

    size_t process(const std::vector<float>& input, std::vector<float>& output) {
        try {
            if (input.empty()) {
                throw std::runtime_error("Input buffer is empty");
            }

            EngineConfig currentConfig;
            {
                std::lock_guard<std::mutex> lock(configMutex);
                currentConfig = config;
            }

            // Ensure input size matches block size
            if (input.size() != currentConfig.blockSize) {
                throw std::runtime_error("Input size must match block size (expected: " + 
                                       std::to_string(currentConfig.blockSize) + 
                                       ", got: " + std::to_string(input.size()) + ")");
            }

            // Check if STN Decomposition is properly initialized
            if (!stnDecomp.isValid()) {
                throw std::runtime_error("STN Decomposition not initialized");
            }

            // Ensure all intermediate buffers are properly sized
            try {
                const size_t bufferSize = currentConfig.blockSize * 2;
                ensureBufferSize(sineComponent, bufferSize);
                ensureBufferSize(transientComponent, bufferSize);
                ensureBufferSize(noiseComponent, bufferSize);
                ensureBufferSize(processedSines, bufferSize);
                ensureBufferSize(processedTransients, bufferSize);
                ensureBufferSize(processedNoise, bufferSize);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Buffer allocation failed: ") + e.what());
            }

            // Step 1: STN Decomposition
            try {
                stnDecomp.decompose(input, sineComponent, transientComponent, noiseComponent);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("STN Decomposition failed: ") + e.what());
            }

            // Step 2: Process each component
            try {
                processComponents(currentConfig);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Component processing failed: ") + e.what());
            }

            // Step 3: Mix components back together
            const size_t outputSize = calculateOutputSize(input.size(), currentConfig);
            try {
                output.resize(outputSize);
                std::fill(output.begin(), output.end(), 0.0f);

                // Mix with appropriate gains
                const float sineGain = 1.0f;
                const float transientGain = 1.0f;
                const float noiseGain = 0.8f;

                for (size_t i = 0; i < outputSize; ++i) {
                    output[i] = sineGain * processedSines[i] +
                               transientGain * processedTransients[i] +
                               noiseGain * processedNoise[i];
                }

                return outputSize;

            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Output mixing failed: ") + e.what());
            }
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Processing failed: ") + e.what());
        }
        return 0; // Should never reach here
    }

    uint32_t getLatency() const {
        std::lock_guard<std::mutex> lock(configMutex);
        return latency;
    }

private:
    void ensureBufferSize(std::vector<float>& buffer, size_t size) {
        if (buffer.size() != size) {
            buffer.resize(size, 0.0f);
        }
    }

    void processComponents(const EngineConfig& currentConfig) {
        // Process sines with pitch shifting
        processSines(currentConfig);
        
        // Process transients with time alignment
        processTransients(currentConfig);
        
        // Process noise with time stretching and spectral morphing
        processNoise(currentConfig);
    }

    void processSines(const EngineConfig& currentConfig) {
        // Create pitch scaler with current config if needed
        if (!pitchScaler || configChanged) {
            pitchScaler = std::make_unique<PitchScalingProcessor>(currentConfig);
        }

        // Process sinusoidal component
        pitchScaler->process(sineComponent, processedSines);

        // Apply time stretching by resampling
        if (std::abs(currentConfig.stretchFactor - 1.0f) > 0.001f) {
            resampleBuffer(processedSines, currentConfig.stretchFactor);
        }
    }

    void processTransients(const EngineConfig& currentConfig) {
        // Preserve transient timing and sharpness
        processedTransients = transientComponent;  // Start with original

        if (std::abs(currentConfig.stretchFactor - 1.0f) > 0.001f) {
            // Time-scale the transient positions
            std::vector<float> timeScaledTransients;
            timeScaledTransients.resize(size_t(processedTransients.size() * currentConfig.stretchFactor));
            std::fill(timeScaledTransients.begin(), timeScaledTransients.end(), 0.0f);

            // Find and preserve transient peaks
            const float threshold = 0.1f;  // Transient detection threshold
            for (size_t i = 1; i < processedTransients.size() - 1; ++i) {
                if (std::abs(processedTransients[i]) > threshold &&
                    std::abs(processedTransients[i]) > std::abs(processedTransients[i-1]) &&
                    std::abs(processedTransients[i]) > std::abs(processedTransients[i+1])) {
                    // Found a transient peak
                    const size_t newPos = size_t(i * currentConfig.stretchFactor);
                    if (newPos < timeScaledTransients.size()) {
                        timeScaledTransients[newPos] = processedTransients[i];
                    }
                }
            }
            
            processedTransients = std::move(timeScaledTransients);
        }
    }

    void processNoise(const EngineConfig& currentConfig) {
        // Create noise morpher with current config if needed
        if (!noiseMorpher || configChanged) {
            noiseMorpher = std::make_unique<NoiseMorphingProcessor>(currentConfig);
        }

        // Process noise component with spectral morphing
        noiseMorpher->process(noiseComponent, processedNoise);
    }

    void resampleBuffer(std::vector<float>& buffer, float stretchFactor) {
        if (buffer.empty() || stretchFactor <= 0.0f) {
            throw std::runtime_error("Invalid buffer or stretch factor");
        }

        const size_t inputSize = buffer.size();
        const size_t outputSize = size_t(inputSize * stretchFactor);
        if (outputSize == 0) {
            throw std::runtime_error("Resample would result in empty buffer");
        }

        std::vector<float> resampled(outputSize);

        for (size_t i = 0; i < outputSize; ++i) {
            const float inputIdx = i / stretchFactor;
            const size_t idx1 = std::min(size_t(inputIdx), inputSize - 1);
            const size_t idx2 = std::min(idx1 + 1, inputSize - 1);
            const float frac = inputIdx - idx1;

            resampled[i] = buffer[idx1] * (1.0f - frac) + buffer[idx2] * frac;
        }

        buffer = std::move(resampled);
    }

    size_t calculateOutputSize(size_t inputSize, const EngineConfig& cfg) const {
        return size_t(inputSize * cfg.stretchFactor);
    }

    // Configuration (ordered to match initialization)
    EngineConfig config;
    STNDecomposition stnDecomp;
    uint32_t latency;
    mutable std::mutex configMutex;
    bool configChanged{false};
    
    // Processing components
    std::unique_ptr<PitchScalingProcessor> pitchScaler;
    std::unique_ptr<NoiseMorphingProcessor> noiseMorpher;

    // Audio buffers (ordered as used)
    std::vector<float> sineComponent;
    std::vector<float> transientComponent;
    std::vector<float> noiseComponent;
    std::vector<float> processedSines;
    std::vector<float> processedTransients;
    std::vector<float> processedNoise;
};

// Public interface implementation
StretchEngine::StretchEngine(const EngineConfig& config)
    : pImpl(std::make_unique<Impl>(config))
{}

StretchEngine::~StretchEngine() = default;

void StretchEngine::updateConfig(const EngineConfig& config) {
    pImpl->updateConfig(config);
}

size_t StretchEngine::process(const std::vector<float>& input,
                            std::vector<float>& output)
{
    return pImpl->process(input, output);
}

uint32_t StretchEngine::getLatency() const {
    return pImpl->getLatency();
}

} // namespace audio
