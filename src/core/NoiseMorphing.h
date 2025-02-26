#pragma once

#include "../../include/stretch_engine.h"
#include <memory>
#include <Eigen/Dense>
#include <fftw3.h>
#include <random>

namespace audio {

class NoiseMorpher {
public:
    explicit NoiseMorpher(const EngineConfig& config);
    ~NoiseMorpher();
    void processFrame(const std::vector<float>& input,
                     std::vector<float>& output,
                     size_t inputOffset,
                     size_t outputOffset);

private:
    Eigen::VectorXf interpolateSpectrum(const Eigen::VectorXf& logMag) const;
    void synthesizeNoise(const Eigen::VectorXf& logMag,
                        size_t outputOffset,
                        std::vector<float>& output);

    const size_t fftSize;
    const size_t hopSize;
    const float stretchFactor;
    
    std::vector<float> window;
    double* fftIn;
    fftw_complex* fftOut;
    fftw_complex* ifftIn;
    double* ifftOut;
    fftw_plan fftPlan;
    fftw_plan ifftPlan;
    
    std::mt19937 rng;
};

class NoiseMorphingProcessor {
public:
    explicit NoiseMorphingProcessor(const EngineConfig& config);
    void process(const std::vector<float>& input, std::vector<float>& output);

private:
    std::unique_ptr<NoiseMorpher> morpher;
    const size_t hopSize;
    const size_t fftSize;
    const EngineConfig config;
};

} // namespace audio
