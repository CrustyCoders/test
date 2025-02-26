# Building and Running the Audio Processing Library

This guide walks you through building the library and processing audio files.

## Prerequisites Installation

### On Ubuntu/Debian:
```bash
# Install required dependencies
sudo apt-get update
sudo apt-get install -y cmake build-essential libfftw3-dev libeigen3-dev libgtest-dev
```

### On macOS:
```bash
# Using Homebrew
brew install cmake fftw eigen googletest
```

### On Windows:
```bash
# Using vcpkg
vcpkg install fftw3:x64-windows eigen3:x64-windows gtest:x64-windows
```

## Building the Library

1. Clone the repository and create a build directory:
```bash
git clone https://github.com/yourusername/stretch_engine.git
cd stretch_engine
mkdir build
cd build
```

2. Configure and build:
```bash
# Configure
cmake ..

# Build (add -j flag for parallel compilation)
cmake --build .
```

3. Run tests to verify everything works:
```bash
ctest --output-on-failure
```

4. Install the library (optional):
```bash
# On Linux/macOS:
sudo cmake --install .

# On Windows (run as administrator):
cmake --install .
```

## Using the Command-Line Tool

The CLI tool accepts WAV files and can perform time-stretching and pitch-shifting operations.

### Example Usage:

1. Double the duration (time stretch factor = 2.0):
```bash
./stretch_engine_cli input.wav output_stretched.wav -s 2.0
```

2. Pitch shift up one octave:
```bash
./stretch_engine_cli input.wav output_pitched.wav -p 12
```

3. Combine time-stretching and pitch-shifting:
```bash
./stretch_engine_cli input.wav output_combined.wav -s 1.5 -p -5
```

4. Customize processing parameters:
```bash
./stretch_engine_cli input.wav output.wav -s 2.0 -p 7 -b 4096 -h 1024 -f 8192
```

### Full Options:
```
Usage: stretch_engine_cli <input.wav> <output.wav> [options]
Options:
  -s <factor>    Time stretch factor (0.125 to 8.0, default: 1.0)
  -p <semitones> Pitch shift in semitones (-12 to +12, default: 0)
  -b <size>      Block size (default: 2048)
  -h <size>      Hop size (default: 512)
  -f <size>      FFT size (default: 4096)
```

## Example Test Files

The repository includes some test audio files in the examples directory:

1. Piano sample (clean harmonic content):
```bash
./stretch_engine_cli examples/piano.wav output_piano.wav -s 2.0 -p 0
```

2. Drum loop (transient-heavy):
```bash
./stretch_engine_cli examples/drums.wav output_drums.wav -s 1.5 -p 0
```

3. Vocal sample (complex harmonic + noise):
```bash
./stretch_engine_cli examples/vocals.wav output_vocals.wav -s 1.0 -p -5
```

## Integrating into Your Project

If you've installed the library, you can use it in your own CMake projects:

```cmake
# In your CMakeLists.txt
find_package(stretch_engine REQUIRED)
target_link_libraries(your_target PRIVATE stretch_engine)
```

Then in your C++ code:

```cpp
#include <stretch_engine.h>
#include <vector>

int main() {
    audio::EngineConfig config;
    config.sampleRate = 44100;
    config.blockSize = 2048;
    config.hopSize = 512;
    config.fftSize = 4096;
    config.stretchFactor = 2.0f;
    config.pitchShift = 0.0f;

    audio::StretchEngine engine(config);

    // Process your audio data
    std::vector<float> input(config.blockSize);
    std::vector<float> output;
    // Fill input with audio data...
    engine.process(input, output);
}
```

## Performance Considerations

- Block size affects latency and CPU usage. Larger blocks = more latency but better CPU efficiency.
- FFT size affects frequency resolution. Larger FFTs = better low-frequency resolution but more latency.
- Hop size affects overlap between frames. Smaller hops = smoother output but more CPU usage.
- For real-time applications, keep block size ≤ 2048 and FFT size ≤ 4096.
