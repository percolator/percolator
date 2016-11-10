include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

file(GLOB PICKED_PROTEIN_SOURCES PickedProteinCaller.cpp Database.cpp Protein.cpp ProteinPeptideIterator.cpp Peptide.cpp PeptideSrc.cpp PeptideConstraint.cpp ../Option.cpp ../Globals.cpp ../MyException.cpp ../Logger.cpp)
add_library(picked_protein STATIC ${PICKED_PROTEIN_SOURCES})
