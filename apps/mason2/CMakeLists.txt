# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/mason2/apps/mason2/CMakeLists.txt
#
# CMakeLists.txt file for mason2.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_mason2_apps_mason2)
message (STATUS "Configuring sandbox/mason2/apps/mason2")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# Enable C++11 in GCC and/or CLANG.
#if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
#    add_definitions ("-std=c++0x")
#endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# We define a library for the reusable parts of Mason.
add_library (mason_sim
             external_split_merge.h
             external_split_merge.cpp
             fragment_generation.h
             genomic_variants.h
             genomic_variants.cpp
             individual_genome.h
             mason_options.cpp
             mason_options.h
             mason_types.h
             methylation_levels.cpp
             methylation_levels.h
             omp_wrapper.h
             sequencing.h
             simulate_454.cpp
             simulate_base.cpp
             simulate_genome.cpp
             simulate_genome.h
             simulate_illumina.cpp
             simulate_sanger.cpp
             variation_size_tsv.h
             vcf_materialization.cpp
             vcf_materialization.h
             )


# Simulate random genomes.
add_executable (mason_genome mason_genome.cpp)

# Simulate methylation levels for gnome.
add_executable (mason_methylation mason_methylation.cpp)

# Simulate sequencing of fragments.
add_executable (mason_frag_sequencing mason_frag_sequencing.cpp)

# Simulation of structural variants into FASTA sequences. 
add_executable (mason_variator mason_variator.cpp)

# The materializer reads in FASTA and VCF and writes out FASTA.
add_executable (mason_materializer mason_materializer.cpp)

# The end-to-end read simulator includes materialization and fragment sampling.
add_executable (mason_simulator mason_simulator.cpp)

# Simulation splicing given a reference, GFF/GTF file, and possibly VCF for variants.
add_executable (mason_splicing mason_splicing.cpp)

# Tests for some of the library functionality.
add_executable (mason_tests mason_tests.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (mason_genome ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_methylation ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_frag_sequencing ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_variator ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_materializer ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_simulator ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_splicing ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_tests ${SEQAN_LIBRARIES} mason_sim)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (mason2)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install mason2 in ${PREFIX}/bin directory
install (TARGETS mason_genome
                 mason_frag_sequencing
                 mason_variator
                 mason_materializer
                 mason_simulator
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/mason2 for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (mason2)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason2")
  set (CPACK_PACKAGE_NAME "mason2")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "mason2")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (mason2 "mason2")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason2")

