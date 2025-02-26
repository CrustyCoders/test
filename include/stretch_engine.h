#ifndef STRETCH_ENGINE_H
#define STRETCH_ENGINE_H

#include <cstdint>
#include <memory>
#include <vector>
#include <complex>

namespace audio {

/**
 * Configuration parameters for the audio processing engine.
 */
struct EngineConfig {
    uint32_t sampleRate{44100};       // Sample rate in Hz
    uint32_t blockSize{2048};         // Processing block size
    uint32_t hopSize{512};            // Hop size for analysis/synthesis
    uint32_t fftSize{4096};           // FFT size for spectral processing
    float stretchFactor{1.0f};        // Time stretch factor (0.125 to 8.0)
    float pitchShift{0.0f};           // Pitch shift in semitones (-12 to +12)
};

/**
 * Main audio processing engine class.
 * Thread-safe for parameter updates during real-time processing.
 */
class StretchEngine {
public:
    /**
     * Creates a new instance of the stretch engine.
     * @param config Initial configuration parameters
     */
    explicit StretchEngine(const EngineConfig& config);
    ~StretchEngine();

    // Prevent copying to ensure thread safety
    StretchEngine(const StretchEngine&) = delete;
    StretchEngine& operator=(const StretchEngine&) = delete;

    /**
     * Updates the engine configuration.
     * Thread-safe: can be called from any thread.
     * @param config New configuration parameters
     */
    void updateConfig(const EngineConfig& config);

    /**
     * Processes a block of audio samples.
     * @param input Input buffer (must be size of blockSize)
     * @param output Output buffer (will be resized as needed)
     * @return Number of output samples generated
     */
    size_t process(const std::vector<float>& input, std::vector<float>& output);

    /**
     * Returns the current latency in samples.
     */
    uint32_t getLatency() const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;  // PIMPL idiom for ABI stability
};

// Forward declaration of STN decomposition component
class STNDecomposition;

} // namespace audio

#endif // STRETCH_ENGINE_H
