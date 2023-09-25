# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.18.2)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(gclib-populate NONE)

include(ExternalProject)
ExternalProject_Add(gclib-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/gpertea/gclib.git"
                    SOURCE_DIR          "/ccb/salz8-1/avaraby/chess3_rerun_31102021/soft/scripts/tiecov_sjs/_deps/gclib-src"
                    BINARY_DIR          "/ccb/salz8-1/avaraby/chess3_rerun_31102021/soft/scripts/tiecov_sjs/_deps/gclib-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


