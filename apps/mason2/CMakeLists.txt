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

# Update the list of file names below if you add source files to your application.
add_executable (mason_genome mason_genome.cpp simulate_genome.h)
add_executable (mason_fragments
                mason_fragments.cpp
                fragment_generation.h
                external_split_merge.h)
add_executable (mason_sequencing
                sequencing.h
                mason_sequencing.cpp
                simulate_base.cpp
                simulate_sanger.cpp
                simulate_454.cpp
                simulate_illumina.cpp)
add_executable (mason_variator mason_variator.cpp variation_size_tsv.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (mason_genome ${SEQAN_LIBRARIES})
target_link_libraries (mason_fragments ${SEQAN_LIBRARIES})
target_link_libraries (mason_sequencing ${SEQAN_LIBRARIES})
target_link_libraries (mason_variator ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (mason2)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install mason2 in ${PREFIX}/bin directory
install (TARGETS mason_genome mason_fragments
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

#seqan_add_app_test (mason2)

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

