#include "../include/stretch_engine.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

namespace {

// Helper function to generate a test sine wave
std::vector<float> generateSineWave(float frequency, float sampleRate, size_t numSamples) {
    std::vector<float> samples(numSamples);
    const float angleIncrement = 2.0f * M_PI * frequency / sampleRate;
    
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i] = std::sin(angleIncrement * i);
    }
    
    return samples;
}

// Helper function to calculate RMS of a signal
float calculateRMS(const std::vector<float>& samples) {
    float sum = 0.0f;
    for (float sample : samples) {
        sum += sample * sample;
    }
    return std::sqrt(sum / samples.size());
}

class StretchEngineTest : public ::testing::Test {
protected:
    void SetUp() override {
        config.sampleRate = 44100;
        config.blockSize = 2048;
        config.hopSize = 512;
        config.fftSize = 4096;
        config.stretchFactor = 1.0f;
        config.pitchShift = 0.0f;
    }

    audio::EngineConfig config;
};

TEST_F(StretchEngineTest, Initialization) {
    ASSERT_NO_THROW({
        audio::StretchEngine engine(config);
    });
}

TEST_F(StretchEngineTest, BasicProcessing) {
    audio::StretchEngine engine(config);
    
    // Generate 1 second of 440 Hz sine wave
    std::vector<float> input = generateSineWave(440.0f, config.sampleRate, config.sampleRate);
    std::vector<float> output;
    
    ASSERT_NO_THROW({
        // Process in blocks
        for (size_t i = 0; i < input.size(); i += config.blockSize) {
            std::vector<float> block(config.blockSize);
            size_t remaining = std::min(static_cast<size_t>(config.blockSize), input.size() - i);
            std::copy(input.begin() + i, input.begin() + i + remaining, block.begin());
            
            if (remaining < config.blockSize) {
                std::fill(block.begin() + remaining, block.end(), 0.0f);
            }
            
            std::vector<float> processedBlock;
            engine.process(block, processedBlock);
            output.insert(output.end(), processedBlock.begin(), processedBlock.end());
        }
    });
    
    // Check that output is non-empty and has reasonable amplitude
    ASSERT_FALSE(output.empty());
    float rms = calculateRMS(output);
    EXPECT_GT(rms, 0.0f);
    EXPECT_LT(rms, 2.0f);
}

TEST_F(StretchEngineTest, TimeStretching) {
    audio::StretchEngine engine(config);
    
    // Generate test signal
    std::vector<float> input = generateSineWave(440.0f, config.sampleRate, config.sampleRate);
    std::vector<float> output;
    
    // Set stretch factor to 2.0 (double duration)
    config.stretchFactor = 2.0f;
    engine.updateConfig(config);
    
    // Process input
    std::vector<float> block(config.blockSize);
    std::copy(input.begin(), input.begin() + config.blockSize, block.begin());
    engine.process(block, output);
    
    // Check output length is approximately double
    EXPECT_NEAR(output.size(), block.size() * 2, config.hopSize);
}

TEST_F(StretchEngineTest, PitchShifting) {
    audio::StretchEngine engine(config);
    
    // Generate test signal
    std::vector<float> input = generateSineWave(440.0f, config.sampleRate, config.sampleRate);
    std::vector<float> output;
    
    // Shift pitch up by 12 semitones (one octave)
    config.pitchShift = 12.0f;
    engine.updateConfig(config);
    
    // Process input
    std::vector<float> block(config.blockSize);
    std::copy(input.begin(), input.begin() + config.blockSize, block.begin());
    engine.process(block, output);
    
    // Check output is non-empty and has reasonable amplitude
    ASSERT_FALSE(output.empty());
    float rms = calculateRMS(output);
    EXPECT_GT(rms, 0.0f);
    EXPECT_LT(rms, 2.0f);
}

TEST_F(StretchEngineTest, CombinedOperations) {
    audio::StretchEngine engine(config);
    
    // Generate test signal
    std::vector<float> input = generateSineWave(440.0f, config.sampleRate, config.sampleRate);
    std::vector<float> output;
    
    // Set both time stretch and pitch shift
    config.stretchFactor = 1.5f;
    config.pitchShift = 7.0f;  // Perfect fifth up
    engine.updateConfig(config);
    
    // Process input
    std::vector<float> block(config.blockSize);
    std::copy(input.begin(), input.begin() + config.blockSize, block.begin());
    engine.process(block, output);
    
    // Check output length and amplitude
    EXPECT_NEAR(output.size(), block.size() * 1.5f, config.hopSize);
    float rms = calculateRMS(output);
    EXPECT_GT(rms, 0.0f);
    EXPECT_LT(rms, 2.0f);
}

TEST_F(StretchEngineTest, Latency) {
    audio::StretchEngine engine(config);
    
    // Check that latency is reported correctly
    uint32_t latency = engine.getLatency();
    EXPECT_EQ(latency, config.fftSize);
    
    // Change FFT size and verify latency updates
    config.fftSize = 8192;
    engine.updateConfig(config);
    latency = engine.getLatency();
    EXPECT_EQ(latency, config.fftSize);
}

} // namespace
