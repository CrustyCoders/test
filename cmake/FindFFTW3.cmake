# Find FFTW3 library
#
# This module defines
# FFTW3_INCLUDE_DIRS, where to find fftw3.h
# FFTW3_LIBRARIES, the libraries to link against
# FFTW3_FOUND, if false, do not try to use FFTW3

find_path(FFTW3_INCLUDE_DIRS NAMES fftw3.h
    PATHS
    /usr/include
    /usr/local/include
)

find_library(FFTW3_LIBRARIES NAMES fftw3
    PATHS
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG
    FFTW3_LIBRARIES
    FFTW3_INCLUDE_DIRS
)

mark_as_advanced(
    FFTW3_INCLUDE_DIRS
    FFTW3_LIBRARIES
)
