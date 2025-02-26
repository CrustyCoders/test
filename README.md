# High-Performance C++ Audio Processing Library

A real-time capable C++ library for high-quality audio time-stretching and pitch-scaling, using advanced Noise Morphing (NM) and Sines-Transients-Noise (STN) decomposition techniques.

## Features

- Modern C++17 implementation
- Real-time processing optimized for DAWs and VST/AU plugins
- High-quality time-stretching (0.125x to 8x) using Noise Morphing
- Artifact-free pitch-shifting (±12 semitones) with formant preservation
- STN decomposition for superior audio quality
- Thread-safe parameter updates
- SIMD optimizations for maximum performance
- Zero allocations on audio thread
- Cross-platform support (Linux, macOS, Windows)

## Requirements

- C++17 compatible compiler
- CMake 3.15 or higher
- FFTW3 library
- Eigen 3.3 or higher
- GoogleTest (for building tests)

### Ubuntu/Debian

```bash
sudo apt-get install cmake build-essential libfftw3-dev libeigen3-dev libgtest-dev
```

### macOS

```bash
brew install cmake fftw eigen googletest
```

### Windows

Using vcpkg:

```bash
vcpkg install fftw3:x64-windows eigen3:x64-windows gtest:x64-windows
```

## Building

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake ..

# Build
cmake --build .

# Run tests
ctest --output-on-failure

# Install
sudo cmake --install .
```

## Usage

### Basic Example

```cpp
#include <stretch_engine.h>
#include <vector>

int main() {
    // Configure the engine
    audio::EngineConfig config;
    config.sampleRate = 44100;
    config.blockSize = 2048;
    config.hopSize = 512;
    config.fftSize = 4096;
    config.stretchFactor = 2.0f;    // Double the duration
    config.pitchShift = -12.0f;     // One octave down

    // Create stretch engine
    audio::StretchEngine engine(config);

    // Process audio blocks
    std::vector<float> input(config.blockSize);
    std::vector<float> output;
    
    // Fill input buffer with audio data...
    engine.process(input, output);
}
```

### CMake Integration

```cmake
find_package(stretch_engine REQUIRED)
target_link_libraries(your_target PRIVATE stretch_engine)
```

## Command-Line Interface

The library includes a command-line tool for processing WAV files:

```bash
stretch_engine_cli input.wav output.wav -s 2.0 -p -12
```

Options:
- `-s <factor>`: Time stretch factor (0.125 to 8.0, default: 1.0)
- `-p <semitones>`: Pitch shift in semitones (-12 to +12, default: 0)
- `-b <size>`: Block size (default: 2048)
- `-h <size>`: Hop size (default: 512)
- `-f <size>`: FFT size (default: 4096)

## Implementation Details

### STN Decomposition

The library uses Sines-Transients-Noise (STN) decomposition to separate the input audio into three components:

1. **Sines**: Harmonic components extracted through spectral peak detection
2. **Transients**: Attack sounds and percussive elements identified by phase deviation
3. **Noise**: Residual energy and textural information

Each component is processed independently with techniques optimized for its characteristics:

- Sinusoidal components use phase vocoder processing with enhanced transient handling
- Transients are preserved and repositioned accurately in time
- Noise is time-stretched using spectral morphing rather than simple phase randomization

### Time-Stretching

The Noise Morphing (NM) algorithm provides high-quality time-stretching by:

- Interpolating log-magnitude spectra for smooth transitions
- Generating white noise excitation signals modulated by the original noise spectra
- Preserving the natural timbre of the sound even at extreme stretch factors

### Pitch-Shifting

Pitch-shifting is achieved through:

- High-quality phase vocoder processing
- Enhanced transient detection and preservation
- Formant preservation for natural-sounding vocal processing

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b new-feature`
3. Make your changes
4. Run the tests: `ctest --output-on-failure`
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.
