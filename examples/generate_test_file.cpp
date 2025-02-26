#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <random>

// Simple WAV header structure (same as in stretch_engine.cpp)
struct WavHeader {
    char riff[4];          // "RIFF"
    uint32_t fileSize;     // File size - 8
    char wave[4];          // "WAVE"
    char fmt[4];           // "fmt "
    uint32_t fmtSize;      // Format chunk size
    uint16_t audioFormat;  // Audio format (1 = PCM)
    uint16_t numChannels;  // Number of channels
    uint32_t sampleRate;   // Sample rate
    uint32_t byteRate;     // Bytes per second
    uint16_t blockAlign;   // Block alignment
    uint16_t bitsPerSample;// Bits per sample
    char data[4];         // "data"
    uint32_t dataSize;    // Data chunk size
};

// Generate a sine wave with specified frequency
std::vector<float> generateSineWave(float frequency, float sampleRate, size_t numSamples) {
    std::vector<float> samples(numSamples);
    float phase = 0.0f;
    float phaseIncrement = 2.0f * M_PI * frequency / sampleRate;
    
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i] = std::sin(phase);
        phase += phaseIncrement;
        if (phase > 2.0f * M_PI) {
            phase -= 2.0f * M_PI;
        }
    }
    
    return samples;
}

// Generate a transient pulse
std::vector<float> generateTransient(size_t numSamples, size_t position) {
    std::vector<float> samples(numSamples, 0.0f);
    
    // Create a sharp attack followed by exponential decay
    if (position < numSamples) {
        samples[position] = 1.0f;
        float decay = 0.95f;
        float amplitude = 1.0f;
        
        for (size_t i = position + 1; i < numSamples && amplitude > 0.001f; ++i) {
            amplitude *= decay;
            samples[i] = amplitude;
        }
    }
    
    return samples;
}

// Generate filtered noise
std::vector<float> generateNoise(size_t numSamples) {
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<float> dist(0.0f, 0.3f);
    
    std::vector<float> samples(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i] = dist(rng);
    }
    
    // Apply simple lowpass filter
    float prevSample = 0.0f;
    const float alpha = 0.2f;
    for (size_t i = 0; i < numSamples; ++i) {
        samples[i] = alpha * samples[i] + (1.0f - alpha) * prevSample;
        prevSample = samples[i];
    }
    
    return samples;
}

void writeWavFile(const std::string& filename,
                 const std::vector<float>& samples,
                 uint32_t sampleRate)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not create file " << filename << std::endl;
        return;
    }

    // Convert float samples to 16-bit
    std::vector<int16_t> rawSamples(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        float sample = samples[i];
        sample = std::max(-1.0f, std::min(1.0f, sample));
        rawSamples[i] = static_cast<int16_t>(sample * 32767.0f);
    }

    WavHeader header;
    std::memcpy(header.riff, "RIFF", 4);
    std::memcpy(header.wave, "WAVE", 4);
    std::memcpy(header.fmt, "fmt ", 4);
    std::memcpy(header.data, "data", 4);

    header.fileSize = static_cast<uint32_t>(sizeof(WavHeader) + rawSamples.size() * 2 - 8);
    header.fmtSize = 16;
    header.audioFormat = 1;
    header.numChannels = 1;
    header.sampleRate = sampleRate;
    header.bitsPerSample = 16;
    header.blockAlign = header.numChannels * header.bitsPerSample / 8;
    header.byteRate = header.sampleRate * header.blockAlign;
    header.dataSize = static_cast<uint32_t>(rawSamples.size() * 2);

    file.write(reinterpret_cast<const char*>(&header), sizeof(WavHeader));
    file.write(reinterpret_cast<const char*>(rawSamples.data()), rawSamples.size() * 2);
}

int main() {
    const uint32_t sampleRate = 44100;
    const size_t duration = 5;  // seconds
    const size_t numSamples = sampleRate * duration;

    // Generate piano-like sample (harmonic content)
    std::vector<float> piano = generateSineWave(440.0f, sampleRate, numSamples);  // A4 note
    // Add harmonics
    for (size_t i = 2; i <= 4; ++i) {
        auto harmonic = generateSineWave(440.0f * i, sampleRate, numSamples);
        for (size_t j = 0; j < numSamples; ++j) {
            piano[j] += harmonic[j] * (0.5f / i);
        }
    }
    writeWavFile("examples/piano.wav", piano, sampleRate);

    // Generate drum-like sample (transients)
    std::vector<float> drums(numSamples, 0.0f);
    for (size_t i = 0; i < duration; ++i) {
        auto hit = generateTransient(numSamples - i * sampleRate, i * sampleRate);
        for (size_t j = 0; j < numSamples - i * sampleRate; ++j) {
            drums[i * sampleRate + j] += hit[j];
        }
    }
    writeWavFile("examples/drums.wav", drums, sampleRate);

    // Generate vocal-like sample (complex harmonic + noise)
    auto fundamental = generateSineWave(220.0f, sampleRate, numSamples);  // Fundamental
    auto formant1 = generateSineWave(600.0f, sampleRate, numSamples);   // First formant
    auto formant2 = generateSineWave(1500.0f, sampleRate, numSamples);  // Second formant
    auto noise = generateNoise(numSamples);  // Breath noise
    
    std::vector<float> vocals(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        vocals[i] = 0.5f * fundamental[i] +
                   0.25f * formant1[i] +
                   0.15f * formant2[i] +
                   0.1f * noise[i];
    }
    writeWavFile("examples/vocals.wav", vocals, sampleRate);

    std::cout << "Generated test files:\n"
              << "  examples/piano.wav\n"
              << "  examples/drums.wav\n"
              << "  examples/vocals.wav\n";

    return 0;
}
