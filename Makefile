# Makefile automates compiling so that OpenCL and MPI are always included and linked

# PATH must contain
# C:/Program Files (x86)/Microsoft Visual Studio/2017/Enterprise/VC/Tools/MSVC/14.11.25503/bin/Hostx86/x86
# so that cl, link and nmake are available on the console.

# Use nmake with this Makefile

# Currently building 32 bit because MPI dll is only 32 bit

# Set the following MSVC, MPI and OpenCL paths to your setup:

# Windows runtime (see [4], [5] and [6] for automating this)
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


################################################################################
# Windows default build
################################################################################

# Link objects
clustermcml-windows.exe: main-windows.o runMCML-windows.o clusterlib-windows.o randomlib-windows.o
	$(MSVC)/link main-windows.o runMCML-windows.o clusterlib-windows.o randomlib-windows.o \
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
main.preprocessed.cpp: main.cpp clusterlib.h
	$(MSVC)/cl main.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"main.preprocessed.cpp"

# Compile the code that should setup and run the CL kernel
runMCML-windows.o: runMCML.preprocessed.cpp
	$(MSVC)/cl runMCML.preprocessed.cpp /c /Fo"runMCML-windows.o"

runMCML.preprocessed.cpp: runMCML.cpp CUDAMCMLio.c randomlib.h Boundary.h Layer.h PhotonTracker.h
	$(MSVC)/cl runMCML.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"runMCML.preprocessed.cpp"

# Compile own libs
clusterlib-windows.o: clusterlib.preprocessed.cpp
	$(MSVC)/cl clusterlib.preprocessed.cpp /c /Fo"clusterlib-windows.o"

clusterlib.preprocessed.cpp: clusterlib.cpp
	$(MSVC)/cl clusterlib.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"clusterlib.preprocessed.cpp"

randomlib-windows.o: randomlib.preprocessed.c
	$(MSVC)/cl randomlib.preprocessed.c /c /Fo"randomlib-windows.o"

randomlib.preprocessed.c: randomlib.c
	$(MSVC)/cl randomlib.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/FI"stdint.h" \
		/P /Fi"randomlib.preprocessed.c"



################################################################################
# Windows CPU build, with debug information
################################################################################

cpumcml-windows.exe: cpumcml-main-windows.o cpumcml-runMCML-windows.o cpumcml-kernel-windows.o randomlib-windows.o
	$(MSVC)/link cpumcml-main-windows.o cpumcml-runMCML-windows.o cpumcml-kernel-windows.o randomlib-windows.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"cpumcml-windows.exe"

cpumcml-main-windows.o: main.cpp clusterlib.h
	$(MSVC)/cl main.cpp /c /Zi \
		/D"CL2CPU" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cpumcml-main-windows.o"

cpumcml-runMCML-windows.o: runMCML.cpp CUDAMCMLio.c randomlib.h Boundary.h Layer.h PhotonTracker.h
	$(MSVC)/cl runMCML.cpp /c /Zi \
		/D"CL2CPU" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cpumcml-runMCML-windows.o"

cpumcml-kernel-windows.o: kernel.c.preprocessed.cpp
	$(MSVC)/cl kernel.c.preprocessed.cpp /c /Zi /Fo"cpumcml-kernel-windows.o"

kernel.c.preprocessed.cpp: kernel.c.cpp Boundary.h Layer.h PhotonTracker.h randomlib.h
	$(MSVC)/cl kernel.c.cpp /c \
		/D"CL2CPU" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/P /Fi"kernel.c.preprocessed.cpp"

kernel.c.cpp: kernel.c cl2cpp.exe
	cl2cpp kernel.c

cl2cpp.exe: cl2cpp.o
	$(MSVC)/link cl2cpp.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"cl2cpp.exe"

cl2cpp.o: cl2cpp.cpp
	$(MSVC)/cl cl2cpp.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cl2cpp.o"


################################################################################
# Windows build, with debug information
# (e.g for [12] Visual Studio, Dr. Memory, AMD CodeXL, Nvidia NSight (?), ...)
################################################################################

# For debug build, pass this rule explicitly to make
clustermcml-windows-debug.exe: main-windows-debug.o runMCML-windows-debug.o clusterlib-windows-debug.o randomlib-windows-debug.o
	$(MSVC)/link main-windows-debug.o runMCML-windows-debug.o clusterlib-windows-debug.o randomlib-windows-debug.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcml-windows-debug.exe"
# Compile new object and generate debug database in separate file
# (see [8] and [9] for explanation of debugging information)
main-windows-debug.o: main.preprocessed.cpp
	$(MSVC)/cl main.preprocessed.cpp /c /Zi /Fo"main-windows-debug.o"
runMCML-windows-debug.o: runMCML.preprocessed.cpp
	$(MSVC)/cl runMCML.preprocessed.cpp /c /Zi /Fo"runMCML-windows-debug.o"
clusterlib-windows-debug.o: clusterlib.preprocessed.cpp
	$(MSVC)/cl clusterlib.preprocessed.cpp /c /Zi /Fo"clusterlib-windows-debug.o"
randomlib-windows-debug.o: randomlib.preprocessed.c
	$(MSVC)/cl randomlib.preprocessed.c /c /Zi /Fo"randomlib-windows-debug.o"



################################################################################
# Windows OpenGL build (for instant visualization of GPU buffers)
# (note that this makes only sense for single GPU)
################################################################################

clustermcml-gl-windows.exe: main-gl-windows.o clusterlib-gl-windows.o randomlib-windows.o runMCML-windows.o gl-windows.o glad.o
	$(MSVC)/link main-gl-windows.o runMCML-windows.o clusterlib-gl-windows.o randomlib-windows.o gl-windows.o glad.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) "user32.lib" "gdi32.lib" "opengl32.lib" \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcml-gl-windows.exe"

main-gl-windows.o: main-gl.preprocessed.cpp
	$(MSVC)/cl main-gl.preprocessed.cpp /c /Fo"main-gl-windows.o"
main-gl.preprocessed.cpp: main.cpp clusterlib.h
	$(MSVC)/cl main.cpp /c \
		/D"GL_VISUALIZATION" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"main-gl.preprocessed.cpp"

clusterlib-gl-windows.o: clusterlib-gl.preprocessed.cpp
	$(MSVC)/cl clusterlib-gl.preprocessed.cpp /c /Fo"clusterlib-gl-windows.o"
clusterlib-gl.preprocessed.cpp: clusterlib.cpp
	$(MSVC)/cl clusterlib.cpp /c \
		/D"GL_VISUALIZATION" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"clusterlib-gl.preprocessed.cpp"


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


################################################################################
# Windows build of Monte Carlo Pi example
################################################################################

# Link objects
clustermcpi-windows.exe: main-windows.o runMonteCarloPi-windows.o clusterlib-windows.o
	$(MSVC)/link main-windows.o runMonteCarloPi-windows.o clusterlib-windows.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustermcpi-windows.exe"

# Compile the code that should setup and run the CL kernel
runMonteCarloPi-windows.o: runMonteCarloPi.preprocessed.cpp
	$(MSVC)/cl runMonteCarloPi.preprocessed.cpp /c /Fo"runMonteCarloPi-windows.o"

runMonteCarloPi.preprocessed.cpp: runMonteCarloPi.cpp
	$(MSVC)/cl runMonteCarloPi.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"runMonteCarloPi.preprocessed.cpp"


################################################################################
# Windows build of Simpson example
################################################################################

# Link objects
clustersimpson-windows.exe: main-windows.o runSimpson-windows.o clusterlib-windows.o
	$(MSVC)/link main-windows.o runSimpson-windows.o clusterlib-windows.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:$(CL_LIBDIR) $(CL_LIBFILE) \
		/LIBPATH:$(MPI_LIBDIR) $(MPI_LIBFILE) \
		/OUT:"clustersimpson-windows.exe"

# Compile the code that should setup and run the CL kernel
runSimpson-windows.o: runSimpson.preprocessed.cpp
	$(MSVC)/cl runSimpson.preprocessed.cpp /c /Fo"runSimpson-windows.o"

runSimpson.preprocessed.cpp: runSimpson.cpp
	$(MSVC)/cl runSimpson.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I$(CL_INCLUDE) /FI$(CL_HEADER) \
		/I$(MPI_INCLUDE) /FI$(MPI_HEADER) \
		/P /Fi"runSimpson.preprocessed.cpp"



clean:
	del *.exe *.o *.preprocessed.c* *.c.cpp *.bin.* *.pdb *.ilk



################################################################################
# Interesting stuff about MSVC:
################################################################################

# [1] https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-by-category

# [2] https://docs.microsoft.com/en-us/cpp/build/reference/linker-options

# [3] https://aras-p.info/blog/2017/10/23/Best-unknown-MSVC-flag-d2cgsummary/

# Automatically find lib and include paths for MSVC runtime:
# [4] https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx#Anchor_1
# [5] https://gist.github.com/FelixK15/3be6a2779c8d3f1f1c354d480fb5cc61
# [6] https://gist.github.com/Kalinovcic/b4d9cc55a37f929cb62320763e8fbb47

# [7] https://lefticus.gitbooks.io/cpp-best-practices/content/02-Use_the_Tools_Available.html#compilers

# [8] https://stackoverflow.com/questions/4659754/the-gs-g-option-equivalent-to-vs2010-cl-compiler

# [9] https://zeuxcg.org/2010/11/22/z7-everything-old-is-new-again/

# [10] https://randomascii.wordpress.com/2014/03/22/make-vc-compiles-fast-through-parallel-compilation/

# [11] https://docs.microsoft.com/en-us/windows/desktop/debug/debug-help-library
# C:\Program Files (x86)\Windows Kits\10\Debuggers

# [12] https://docs.microsoft.com/en-us/visualstudio/debugger/how-to-debug-an-executable-not-part-of-a-visual-studio-solution?view=vs-2017
