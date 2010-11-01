link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB ELUDE_SOURCES *.cpp)

add_library(eludelibrary SHARED RetentionFeatures.cpp DataManager.cpp EludeMain.cpp LibSVRModel.cpp LibsvmWrapper.cpp SVRModel.h RetentionModel.cpp EludeCaller.cpp  LTSRegression.cpp ../svm.cpp ../Normalizer.cpp ../UniNormalizer.cpp ../StdvNormalizer.cpp ../Option.cpp ../Enzyme.cpp ../PSMDescription.cpp ../Globals.cpp)

add_executable(elude ${ELUDE_SOURCES})

target_link_libraries(elude perclibrary ArrayLibrary  ${XERCESC_LIBRARIES} )

install(TARGETS elude EXPORT PERCOLATOR DESTINATION .) # Important to use relative path here (used by CPack)!

IF(STATIC)
  set_property(TARGET elude PROPERTY LINK_SEARCH_END_STATIC ON)
ENDIF()
