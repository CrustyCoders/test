#pragma once

#include "../../include/stretch_engine.h"
#include <memory>
#include <vector>

namespace audio {

// Forward declaration for implementation details
namespace detail {
    class STNDecompositionImpl;
}

class STNDecomposition {
public:
    explicit STNDecomposition(const EngineConfig& config);
    ~STNDecomposition();

    // Deleted copy operations but allow moves
    STNDecomposition(const STNDecomposition&) = delete;
    STNDecomposition& operator=(const STNDecomposition&) = delete;
    STNDecomposition(STNDecomposition&&) noexcept = default;
    STNDecomposition& operator=(STNDecomposition&&) noexcept = default;

    void decompose(const std::vector<float>& input,
                  std::vector<float>& sines,
                  std::vector<float>& transients,
                  std::vector<float>& noise);

    // Check if the implementation is properly initialized
    bool isValid() const { return pImpl != nullptr; }

private:
    std::unique_ptr<detail::STNDecompositionImpl> pImpl;
};

} // namespace audio
