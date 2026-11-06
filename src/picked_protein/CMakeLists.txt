include_directories(${PERCOLATOR_SOURCE_DIR}/src)

add_library(picked_protein STATIC PickedProteinCaller.cpp Database.cpp Protein.cpp ProteinPeptideIterator.cpp Peptide.cpp PeptideSrc.cpp PeptideConstraint.cpp ProteinDigest.cpp)
# Globals (verbosity) comes from perclibrary, which in turn uses this library
# through PickedProteinInterface.
target_link_libraries(picked_protein PUBLIC perclibrary)
