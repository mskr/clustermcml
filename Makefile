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

mcml-comparison.exe: mcmlio.o mcmlnr.o mcmlgo.o main.o cpumcml-runMCML-windows.o cpumcml-kernel-windows.o
	$(MSVC)/link main.o mcmlio.o mcmlnr.o mcmlgo.o cpumcml-runMCML-windows.o cpumcml-kernel-windows.o \
		/DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"mcml-comparison.exe"

main.o: main.c
	$(MSVC)/cl main.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"main.o"

mcmlio.o: mcmlio.c
	$(MSVC)/cl mcmlio.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"mcmlio.o"

mcmlnr.o: mcmlnr.c
	$(MSVC)/cl mcmlnr.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"mcmlnr.o"

mcmlgo.o: mcmlgo.cpp
	$(MSVC)/cl mcmlgo.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"mcmlgo.o"



cpu-main.o: cpu-main.cpp
	$(MSVC)/cl cpu-main.cpp /c /Zi \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cpu-main.o"



cpumcml-runMCML-windows.o: runMCML.cpp
	$(MSVC)/cl runMCML.cpp /c /Zi \
		/D"CL2CPU" \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cpumcml-runMCML-windows.o"


cpumcml-kernel-windows.o: kernel.c.preprocessed.cpp
	$(MSVC)/cl kernel.c.preprocessed.cpp /c /Zi /Fo"cpumcml-kernel-windows.o"

kernel.c.preprocessed.cpp: kernel.c.cpp
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




clean:
	rm *.preprocessed.cpp *.exe *.o *.c.cpp
