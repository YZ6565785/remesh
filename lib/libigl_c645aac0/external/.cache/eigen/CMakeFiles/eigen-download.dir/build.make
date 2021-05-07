# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen

# Utility rule file for eigen-download.

# Include the progress variables for this target.
include CMakeFiles/eigen-download.dir/progress.make

CMakeFiles/eigen-download: CMakeFiles/eigen-download-complete


CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-install
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-mkdir
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-download
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-update
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-patch
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-build
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-install
CMakeFiles/eigen-download-complete: eigen-download-prefix/src/eigen-download-stamp/eigen-download-test
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'eigen-download'"
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles
	/usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles/eigen-download-complete
	/usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-done

eigen-download-prefix/src/eigen-download-stamp/eigen-download-install: eigen-download-prefix/src/eigen-download-stamp/eigen-download-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "No install step for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E echo_append
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-install

eigen-download-prefix/src/eigen-download-stamp/eigen-download-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Creating directories for 'eigen-download'"
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/cmake/../external/eigen
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/tmp
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp
	/usr/bin/cmake -E make_directory /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src
	/usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-mkdir

eigen-download-prefix/src/eigen-download-stamp/eigen-download-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-gitinfo.txt
eigen-download-prefix/src/eigen-download-stamp/eigen-download-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (git clone) for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external && /usr/bin/cmake -P /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/tmp/eigen-download-gitclone.cmake
	cd /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external && /usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-download

eigen-download-prefix/src/eigen-download-stamp/eigen-download-update: eigen-download-prefix/src/eigen-download-stamp/eigen-download-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Performing update step for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/eigen && /usr/bin/cmake -P /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/tmp/eigen-download-gitupdate.cmake

eigen-download-prefix/src/eigen-download-stamp/eigen-download-patch: eigen-download-prefix/src/eigen-download-stamp/eigen-download-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "No patch step for 'eigen-download'"
	/usr/bin/cmake -E echo_append
	/usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-patch

eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure: eigen-download-prefix/tmp/eigen-download-cfgcmd.txt
eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure: eigen-download-prefix/src/eigen-download-stamp/eigen-download-update
eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure: eigen-download-prefix/src/eigen-download-stamp/eigen-download-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No configure step for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E echo_append
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure

eigen-download-prefix/src/eigen-download-stamp/eigen-download-build: eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No build step for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E echo_append
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-build

eigen-download-prefix/src/eigen-download-stamp/eigen-download-test: eigen-download-prefix/src/eigen-download-stamp/eigen-download-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "No test step for 'eigen-download'"
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E echo_append
	cd /home/yuhang/Desktop/COMP0119/final_project/build/eigen-build && /usr/bin/cmake -E touch /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/eigen-download-prefix/src/eigen-download-stamp/eigen-download-test

eigen-download: CMakeFiles/eigen-download
eigen-download: CMakeFiles/eigen-download-complete
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-install
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-mkdir
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-download
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-update
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-patch
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-configure
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-build
eigen-download: eigen-download-prefix/src/eigen-download-stamp/eigen-download-test
eigen-download: CMakeFiles/eigen-download.dir/build.make

.PHONY : eigen-download

# Rule to build all files generated by this target.
CMakeFiles/eigen-download.dir/build: eigen-download

.PHONY : CMakeFiles/eigen-download.dir/build

CMakeFiles/eigen-download.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/eigen-download.dir/cmake_clean.cmake
.PHONY : CMakeFiles/eigen-download.dir/clean

CMakeFiles/eigen-download.dir/depend:
	cd /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen /home/yuhang/Desktop/COMP0119/final_project/lib/libigl_c645aac0/external/.cache/eigen/CMakeFiles/eigen-download.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/eigen-download.dir/depend

