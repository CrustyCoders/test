#pragma once

#include "STNDecompositionDecl.h"
#include <Eigen/Dense>
#include <fftw3.h>

namespace audio {
namespace detail {

    class STNDecompositionImpl {
    public:
        STNDecompositionImpl(const EngineConfig& config);
        ~STNDecompositionImpl();
        
        // Implement move operations
        STNDecompositionImpl(STNDecompositionImpl&& other) noexcept;
        STNDecompositionImpl& operator=(STNDecompositionImpl&& other) noexcept;
        
        // Delete copy operations
        STNDecompositionImpl(const STNDecompositionImpl&) = delete;
        STNDecompositionImpl& operator=(const STNDecompositionImpl&) = delete;
    
    void decompose(const std::vector<float>& input,
                  std::vector<float>& sines,
                  std::vector<float>& transients,
                  std::vector<float>& noise);

private:
    void decomposeSpectrum(const Eigen::VectorXf& magnitude,
                          const Eigen::VectorXf& phase,
                          size_t frameOffset);

    const size_t fftSize;
    const size_t hopSize;
    
    std::vector<float> window;
    double* fftIn;
    fftw_complex* fftOut;
    fftw_complex* ifftIn;
    double* ifftOut;
    fftw_plan fftPlan;
    fftw_plan ifftPlan;

    Eigen::VectorXf sinesMagnitude;
    Eigen::VectorXf sinesPhase;
    Eigen::VectorXf transientsMagnitude;
    Eigen::VectorXf transientsPhase;
    Eigen::VectorXf noiseMagnitude;
    Eigen::VectorXf noisePhase;
};

} // namespace detail
} // namespace audio
