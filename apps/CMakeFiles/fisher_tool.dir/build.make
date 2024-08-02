# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /opt/gw_analysis_tools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /opt/gw_analysis_tools

# Include any dependencies generated for this target.
include apps/CMakeFiles/fisher_tool.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/fisher_tool.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/fisher_tool.dir/flags.make

apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o: apps/CMakeFiles/fisher_tool.dir/flags.make
apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o: apps/fisher_tool.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/opt/gw_analysis_tools/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o"
	cd /opt/gw_analysis_tools/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o -c /opt/gw_analysis_tools/apps/fisher_tool.cpp

apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fisher_tool.dir/fisher_tool.cpp.i"
	cd /opt/gw_analysis_tools/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /opt/gw_analysis_tools/apps/fisher_tool.cpp > CMakeFiles/fisher_tool.dir/fisher_tool.cpp.i

apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fisher_tool.dir/fisher_tool.cpp.s"
	cd /opt/gw_analysis_tools/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /opt/gw_analysis_tools/apps/fisher_tool.cpp -o CMakeFiles/fisher_tool.dir/fisher_tool.cpp.s

# Object files for target fisher_tool
fisher_tool_OBJECTS = \
"CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o"

# External object files for target fisher_tool
fisher_tool_EXTERNAL_OBJECTS =

bin/fisher_tool: apps/CMakeFiles/fisher_tool.dir/fisher_tool.cpp.o
bin/fisher_tool: apps/CMakeFiles/fisher_tool.dir/build.make
bin/fisher_tool: lib/libgwat.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libpthread.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libsz.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libz.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libdl.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libm.so
bin/fisher_tool: /usr/local/lib64/libadolc.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libgsl.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libgslcblas.so
bin/fisher_tool: /usr/local/lib/libbayesship.so
bin/fisher_tool: /usr/lib/gcc/x86_64-linux-gnu/10/libgomp.so
bin/fisher_tool: /usr/lib/x86_64-linux-gnu/libpthread.so
bin/fisher_tool: apps/CMakeFiles/fisher_tool.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/opt/gw_analysis_tools/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/fisher_tool"
	cd /opt/gw_analysis_tools/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fisher_tool.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/fisher_tool.dir/build: bin/fisher_tool

.PHONY : apps/CMakeFiles/fisher_tool.dir/build

apps/CMakeFiles/fisher_tool.dir/clean:
	cd /opt/gw_analysis_tools/apps && $(CMAKE_COMMAND) -P CMakeFiles/fisher_tool.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/fisher_tool.dir/clean

apps/CMakeFiles/fisher_tool.dir/depend:
	cd /opt/gw_analysis_tools && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/gw_analysis_tools /opt/gw_analysis_tools/apps /opt/gw_analysis_tools /opt/gw_analysis_tools/apps /opt/gw_analysis_tools/apps/CMakeFiles/fisher_tool.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/fisher_tool.dir/depend

