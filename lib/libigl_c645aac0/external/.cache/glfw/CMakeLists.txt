# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 2.8.2)

project(glfw-download NONE)

include(ExternalProject)
ExternalProject_Add(glfw-download
                    TLS_VERIFY;OFF;GIT_CONFIG;advice.detachedHead=false;GIT_REPOSITORY;https://github.com/glfw/glfw.git;GIT_TAG;3.3
                    SOURCE_DIR          "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw"
                    BINARY_DIR          "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/build/glfw-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
