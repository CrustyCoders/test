#include "../../include/stretch_engine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

// Simple WAV file header structure
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

bool readWavFile(const std::string& filename,
                std::vector<float>& samples,
                uint32_t& sampleRate)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    WavHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(WavHeader));

    if (std::strncmp(header.riff, "RIFF", 4) != 0 ||
        std::strncmp(header.wave, "WAVE", 4) != 0) {
        std::cerr << "Error: Invalid WAV file format\n";
        return false;
    }

    if (header.audioFormat != 1) {
        std::cerr << "Error: Only PCM format is supported\n";
        return false;
    }

    if (header.bitsPerSample != 16) {
        std::cerr << "Error: Only 16-bit audio is supported\n";
        return false;
    }

    sampleRate = header.sampleRate;
    const size_t numSamples = header.dataSize / (header.bitsPerSample / 8);
    samples.resize(numSamples);

    // Read 16-bit samples and convert to float
    std::vector<int16_t> rawSamples(numSamples);
    file.read(reinterpret_cast<char*>(rawSamples.data()), header.dataSize);

    for (size_t i = 0; i < numSamples; ++i) {
        samples[i] = rawSamples[i] / 32768.0f;
    }

    return true;
}

bool writeWavFile(const std::string& filename,
                 const std::vector<float>& samples,
                 uint32_t sampleRate)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not create file " << filename << std::endl;
        return false;
    }

    // Convert float samples to 16-bit
    std::vector<int16_t> rawSamples(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        float sample = samples[i];
        // Clamp to [-1, 1]
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

    return true;
}

void printUsage() {
    std::cout << "Usage: stretch_engine <input.wav> <output.wav> [options]\n"
              << "Options:\n"
              << "  -s <factor>    Time stretch factor (0.125 to 8.0, default: 1.0)\n"
              << "  -p <semitones> Pitch shift in semitones (-12 to +12, default: 0)\n"
              << "  -b <size>      Block size (default: 2048)\n"
              << "  -h <size>      Hop size (default: 512)\n"
              << "  -f <size>      FFT size (default: 4096)\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage();
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    // Default configuration
    audio::EngineConfig config;
    config.stretchFactor = 1.0f;
    config.pitchShift = 0.0f;

    // Parse command line options
    for (int i = 3; i < argc; i += 2) {
        if (i + 1 >= argc) {
            std::cerr << "Error: Missing value for option " << argv[i] << std::endl;
            return 1;
        }

        std::string option = argv[i];
        float value = std::stof(argv[i + 1]);

        if (option == "-s") {
            config.stretchFactor = std::max(0.125f, std::min(8.0f, value));
        } else if (option == "-p") {
            config.pitchShift = std::max(-12.0f, std::min(12.0f, value));
        } else if (option == "-b") {
            config.blockSize = static_cast<uint32_t>(value);
        } else if (option == "-h") {
            config.hopSize = static_cast<uint32_t>(value);
        } else if (option == "-f") {
            config.fftSize = static_cast<uint32_t>(value);
        } else {
            std::cerr << "Unknown option: " << option << std::endl;
            return 1;
        }
    }

    // Read input file
    std::vector<float> inputSamples;
    uint32_t sampleRate;
    if (!readWavFile(inputFile, inputSamples, sampleRate)) {
        return 1;
    }
    config.sampleRate = sampleRate;

    try {
        std::cout << "Creating stretch engine with config:\n"
                  << "  Sample rate: " << config.sampleRate << " Hz\n"
                  << "  Block size: " << config.blockSize << "\n"
                  << "  Hop size: " << config.hopSize << "\n"
                  << "  FFT size: " << config.fftSize << "\n"
                  << "  Stretch factor: " << config.stretchFactor << "\n"
                  << "  Pitch shift: " << config.pitchShift << " semitones\n";

        // Create and configure the stretch engine
        audio::StretchEngine engine(config);
        std::cout << "Engine created successfully\n";

        // Process the audio in blocks
        std::vector<float> output;
        std::vector<float> block;
        block.resize(config.blockSize);

        size_t processedSamples = 0;
        std::vector<float> processedOutput;
        const size_t totalBlocks = (inputSamples.size() + config.blockSize - 1) / config.blockSize;
        size_t currentBlock = 0;

        std::cout << "Processing " << inputSamples.size() << " samples in " 
                  << totalBlocks << " blocks...\n";

        while (processedSamples < inputSamples.size()) {
            // Prepare input block
            size_t remainingSamples = inputSamples.size() - processedSamples;
            size_t samplesToProcess = std::min(remainingSamples, static_cast<size_t>(config.blockSize));
            
            if (samplesToProcess < config.blockSize) {
                // Zero-pad the last block
                std::fill(block.begin(), block.end(), 0.0f);
            }
            
            std::copy(inputSamples.begin() + processedSamples,
                     inputSamples.begin() + processedSamples + samplesToProcess,
                     block.begin());

            // Process block
            std::cout << "Processing block " << ++currentBlock << "/" << totalBlocks 
                      << " (" << samplesToProcess << " samples)\n";
            engine.process(block, output);

            // Accumulate output
            std::cout << "Block produced " << output.size() << " output samples\n";
            processedOutput.insert(processedOutput.end(),
                                output.begin(),
                                output.end());

            processedSamples += config.blockSize;
        }

        std::cout << "Processing complete. Writing " << processedOutput.size() 
                  << " samples to output file...\n";

        // Write output file
        if (!writeWavFile(outputFile, processedOutput, sampleRate)) {
            std::cerr << "Failed to write output file\n";
            return 1;
        }

        std::cout << "Processing complete. Output written to " << outputFile << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}
