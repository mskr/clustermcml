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

msvc-mcml.exe: msvc-mcmlio.o msvc-mcmlnr.o msvc-mcmlgo.o msvc-mcmlmain.o cpu-runMCML.o cpu-kernel.o
	$(MSVC)/link msvc-mcmlmain.o msvc-mcmlio.o msvc-mcmlnr.o msvc-mcmlgo.o cpu-runMCML.o cpu-kernel.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"msvc-mcml.exe"

msvc-mcmlmain.o: mcmlmain.c
	$(MSVC)/cl mcmlmain.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"msvc-mcmlmain.o"

msvc-mcmlio.o: mcmlio.c
	$(MSVC)/cl mcmlio.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"msvc-mcmlio.o"

msvc-mcmlnr.o: mcmlnr.c
	$(MSVC)/cl mcmlnr.c /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Zi \
		/c /Fo"msvc-mcmlnr.o"

msvc-mcmlgo.o: mcmlgo.cpp
	$(MSVC)/cl mcmlgo.cpp /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \



cpu-main.o: cpu-main.cpp
	$(MSVC)/cl cpu-main.cpp /c /Zi \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cpu-main.o"







clean:
	rm *.preprocessed.cpp *.exe *.o *.c.cpp
