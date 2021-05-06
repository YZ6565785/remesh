
if(NOT "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitinfo.txt" IS_NEWER_THAN "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/opt/homebrew/bin/git" -c http.sslVerify=false clone --no-checkout --config "advice.detachedHead=false" "https://github.com/glfw/glfw.git" "glfw"
    WORKING_DIRECTORY "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/glfw/glfw.git'")
endif()

execute_process(
  COMMAND "/opt/homebrew/bin/git" -c http.sslVerify=false checkout 3.3 --
  WORKING_DIRECTORY "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '3.3'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/opt/homebrew/bin/git" -c http.sslVerify=false submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/cmake/../external/glfw'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitinfo.txt"
    "/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/yuhangzhang/Desktop/UCL/courses/Term2/COMP0119_Acquisition_and_Processing_of_3D_Geometry/final_project/lib/libigl_c645aac0/external/.cache/glfw/glfw-download-prefix/src/glfw-download-stamp/glfw-download-gitclone-lastrun.txt'")
endif()

