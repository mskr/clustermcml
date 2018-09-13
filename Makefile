# Makefile automates compiling so that OpenCL and MPI are always included and linked

# PATH must contain
# C:/Program Files (x86)/Microsoft Visual Studio/2017/Enterprise/VC/Tools/MSVC/14.11.25503/bin/Hostx86/x86
# so that cl, link and nmake are available on the console.

# Use nmake with this Makefile

# Currently building 32 bit because MPI dll is only 32 bit

# Set the following MSVC, MPI and OpenCL paths to your setup:

# Windows runtime
MSVC = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128\bin\Hostx86\x86"
MSVC_INCLUDE = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128\include"
MSVC_INCLUDE_UCRT = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\ucrt"
MSVC_INCLUDE_UM = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\um"
MSVC_INCLUDE_SHARED = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\shared"
MSVC_LIB = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128/lib/x86"
MSVC_LIB_UCRT = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.16299.0\ucrt\x86"
MSVC_LIB_UM = "C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86"

# OpenCL
CL_INCLUDE = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\include"
CL_HEADER = "CL/opencl.h"
CL_LIBDIR = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\Win32"
CL_LIBFILE = "OpenCL.lib"

# MPI
MPI_INCLUDE = "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
MPI_HEADER = "mpi.h"
MPI_LIBDIR = "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86"
MPI_LIBFILE = "msmpi.lib"

cpu-mcml.exe: cpu-mcml.o
	$(MSVC)/link cpu-mcml.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"cpu-mcml.exe"

cpu-mcml.o: cpu-mcml.cpp kernel.c
	$(MSVC)/cl cpu-mcml.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/D"CL2CPP" \
		/c /Fo"cpu-mcml.o"

# Link objects
clustermcml-windows.exe: main-windows.o runMCML-windows.o
	$(MSVC)/link main-windows.o runMCML-windows.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcml-windows.exe"

# Compile new object
main-windows.o: main.preprocessed.cpp
	$(MSVC)/cl main.preprocessed.cpp /c /Fo"main-windows.o"

# Preprocess source
main.preprocessed.cpp: main.cpp
	$(MSVC)/cl main.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"main.preprocessed.cpp"

# Compile the code that should setup and run the CL kernel
runMCML-windows.o: runMCML.preprocessed.cpp
	$(MSVC)/cl runMCML.preprocessed.cpp /c /Fo"runMCML-windows.o"

runMCML.preprocessed.cpp: runMCML.cpp CUDAMCMLio.c
	$(MSVC)/cl runMCML.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"runMCML.preprocessed.cpp"


##################################################################################################
# Windows, with debug information (e.g for Dr. Memory, AMD CodeXL, Nvidia NSight (?), ...)
##################################################################################################

# For debug build, pass this rule explicitly to make
clustermcml-windows-debug.exe: main-windows-debug.o runMCML-windows-debug.o
	$(MSVC)/link main-windows-debug.o runMCML-windows-debug.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcml-windows-debug.exe"
# Compile new object and generate debug database in separate file
main-windows-debug.o: main.preprocessed.cpp
	$(MSVC)/cl main.preprocessed.cpp /c /Zi /Fo"main-windows-debug.o"
runMCML-windows-debug.o: runMCML.preprocessed.cpp
	$(MSVC)/cl runMCML.preprocessed.cpp /c /Zi /Fo"runMCML-windows-debug.o"



######################################################
# Windows, with OpenGL (for instant visualization)
######################################################

clustermcml-gl-windows.exe: main-gl-windows.o runMCML-windows.o gl-windows.o glad.o
	$(MSVC)/link main-gl-windows.o runMCML-windows.o gl-windows.o glad.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) "user32.lib" "gdi32.lib" "opengl32.lib" \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcml-gl-windows.exe"

main-gl-windows.o: main-gl.preprocessed.cpp
	$(MSVC)/cl main-gl.preprocessed.cpp /c /Fo"main-gl-windows.o"
main-gl.preprocessed.cpp: main.cpp
	$(MSVC)/cl main.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/D"GL_VISUALIZATION" \
		/P /Fi"main-gl.preprocessed.cpp"

gl-windows.o: gl-windows.preprocessed.cpp
	$(MSVC)/cl gl-windows.preprocessed.cpp /c /Fo"gl-windows.o"
gl-windows.preprocessed.cpp: gl-windows.cpp
	$(MSVC)/cl gl-windows.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(MSVC_INCLUDE_UM) \
		/I$(MSVC_INCLUDE_SHARED) \
		/P /Fi"gl-windows.preprocessed.cpp"

# GL loader generated with http://glad.dav1d.de/
glad.o: glad.c
	$(MSVC)/cl glad.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(MSVC_INCLUDE_UM) \
		/I$(MSVC_INCLUDE_SHARED) \
		/Fo"glad.o"


clean:
	del *.exe *.o *.preprocessed.cpp *.cl.bin.* *.pdb *.ilk

# Interesting stuff about MSVC:

# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-by-category

# https://docs.microsoft.com/en-us/cpp/build/reference/linker-options

# https://aras-p.info/blog/2017/10/23/Best-unknown-MSVC-flag-d2cgsummary/

# https://gist.github.com/FelixK15/3be6a2779c8d3f1f1c354d480fb5cc61

# https://lefticus.gitbooks.io/cpp-best-practices/content/02-Use_the_Tools_Available.html#compilers

# https://stackoverflow.com/questions/4659754/the-gs-g-option-equivalent-to-vs2010-cl-compiler

# https://zeuxcg.org/2010/11/22/z7-everything-old-is-new-again/

# https://randomascii.wordpress.com/2014/03/22/make-vc-compiles-fast-through-parallel-compilation/

# https://docs.microsoft.com/en-us/windows/desktop/debug/debug-help-library
# C:\Program Files (x86)\Windows Kits\10\Debuggers