# Makefile automates compiling so that OpenCL and MPI are always included and linked

# PATH must contain
# C:/Program Files (x86)/Microsoft Visual Studio/2017/Enterprise/VC/Tools/MSVC/14.11.25503/bin/Hostx86/x86
# so that cl, link and nmake are available on the console.

# Set MSVC, MPI and OpenCL paths to your setup!

# Use nmake with this Makefile

# Currently building 32 bit because MPI dll is only 32 bit

# Runtime
MSVC_INCLUDE = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\include"
MSVC_INCLUDE_UCRT = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.10240.0\ucrt"
MSVC_LIB = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428/lib/x86"
MSVC_LIB_UCRT = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.10240.0\ucrt\x86"
MSVC_LIB_UM = "C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86"

# Link objects
clustermcml-windows.exe: main-windows.o runMCML-windows.o
	link main-windows.o runMCML-windows.o \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/LIBPATH:"C:\common-lib-amd-APPSDK-3.0-master\3-0\lib\x86" OpenCL.lib \
		/LIBPATH:"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86" msmpi.lib \
		/OUT:"clustermcml-windows.exe"


# Compile new object when source is newer
main-windows.o: main.preprocessed.cpp
	cl main.preprocessed.cpp /c /Fo"main-windows.o"

# Preprocess source
main.preprocessed.cpp: main.cpp
	cl main.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"C:\common-lib-amd-APPSDK-3.0-master\3-0\include" /FI"CL/opencl.h" \
		/I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" /FI"mpi.h" \
		/P /Fi"main.preprocessed.cpp"

# Compile the code that should setup and run the CL kernel
runMCML-windows.o: runMCML.preprocessed.cpp
	cl runMCML.preprocessed.cpp /c /Fo"runMCML-windows.o"

runMCML.preprocessed.cpp: runMCML.cpp CUDAMCMLio.c
	cl runMCML.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"C:\common-lib-amd-APPSDK-3.0-master\3-0\include" /FI"CL/opencl.h" \
		/I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" /FI"mpi.h" \
		/P /Fi"runMCML.preprocessed.cpp"

clean:
	del *-windows.exe *-windows.o *.preprocessed.cpp *.cl.bin.*


# Interesting stuff about MSVC:

# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-by-category

# https://docs.microsoft.com/en-us/cpp/build/reference/linker-options

# https://aras-p.info/blog/2017/10/23/Best-unknown-MSVC-flag-d2cgsummary/

# https://gist.github.com/FelixK15/3be6a2779c8d3f1f1c354d480fb5cc61

# https://lefticus.gitbooks.io/cpp-best-practices/content/02-Use_the_Tools_Available.html#compilers