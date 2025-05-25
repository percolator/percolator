# cmake/patch_eigen.cmake
file(READ "${EIGEN_FILE}" EIGEN_CONTENTS)

# Replace the line that causes the MSVC '_finite' error
string(REPLACE "using std::_finite;" 
"// Workaround for MSVC '_finite' bug
#if defined(_MSC_VER)
using ::_finite;
#else
using std::isfinite;
#endif" 
PATCHED_CONTENTS "${EIGEN_CONTENTS}")

file(WRITE "${EIGEN_FILE}" "${PATCHED_CONTENTS}")
message(STATUS "Applied patch to Eigen: workaround for MSVC '_finite'")
