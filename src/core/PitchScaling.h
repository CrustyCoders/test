#pragma once

#include "../../include/stretch_engine.h"
#include <memory>
#include <Eigen/Dense>
#include <fftw3.h>

namespace audio {

class PitchScaler {
public:
    explicit PitchScaler(const EngineConfig& config);
    ~PitchScaler();
    void processFrame(const std::vector<float>& input,
                     std::vector<float>& output,
                     size_t inputOffset,
                     size_t outputOffset);

private:
    void preserveFormants(const Eigen::VectorXf& originalMagnitude);
    
    const size_t fftSize;
    const size_t hopSize;
    const float pitchFactor;
    const float sampleRate;

    std::vector<float> window;
    std::vector<float> lastPhase;
    std::vector<float> sumPhase;
    
    double* fftIn;
    fftw_complex* fftOut;
    fftw_complex* ifftIn;
    double* ifftOut;
    fftw_plan fftPlan;
    fftw_plan ifftPlan;
};

class PitchScalingProcessor {
public:
    explicit PitchScalingProcessor(const EngineConfig& config);
    void process(const std::vector<float>& input, std::vector<float>& output);

private:
    std::unique_ptr<PitchScaler> scaler;
    const size_t hopSize;
    const size_t fftSize;
};

} // namespace audio
