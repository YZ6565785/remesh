cmake_minimum_required(VERSION 3.1)

project(remeshing)

set(CMAKE_CXX_STANDARD 11)

# message([STATUS] "CMAKE_CURRENT_SOURCE_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")
message([STATUS] "ENV CGAL_DIR : $ENV{CGAL_DIR}")

# set the directory of cmake macro file 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# set libigl's magic options 
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI   "Use GLFW_IMGUI"  ON )

# search libraries
find_package(LIBIGL REQUIRED ) 
find_package(CGAL REQUIRED)

add_executable(remeshing src/main.cpp)

# link library
target_link_libraries(remeshing PUBLIC
  igl::core 
  igl::opengl 
  igl::opengl_glfw 
  igl::opengl_glfw_imgui 
  CGAL::CGAL
)

