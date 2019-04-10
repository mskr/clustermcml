MSVC                = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\bin\Hostx86\x86"
MSVC_INCLUDE        = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\include" # STL
MSVC_INCLUDE_UCRT   = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\ucrt" # crtdefs.h
MSVC_INCLUDE_UM     = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\um"
MSVC_INCLUDE_SHARED = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\shared"
MSVC_LIB            = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\lib\x86"
MSVC_LIB_UCRT       = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\ucrt\x86"
MSVC_LIB_UM         = "C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86"

setupInputArrays.exe: setupInputArrays.o clmem-windows.o
	$(MSVC)/link setupInputArrays.o clmem-windows.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"setupInputArrays.exe"

setupInputArrays.o: setupInputArrays.cpp ../../CUDAMCMLio.h ../../CUDAMCMLio.c ../../clcheck.h ../../clmem.h ../../geometrylib.h ../../Layer.h ../../Boundary.h
	$(MSVC)/cl setupInputArrays.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/I.. \
		/D"NO_GPU" \
		/c /Fo"setupInputArrays.o"

clmem-windows.o: clmem.preprocessed.cpp
	$(MSVC)/cl clmem.preprocessed.cpp /c /W3 /Fo"clmem-windows.o"

clmem.preprocessed.cpp: ../../clmem.cpp
	$(MSVC)/cl "../../clmem.cpp" /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/D"NO_GPU" \
		/P /Fi"clmem.preprocessed.cpp"


#####################################################################

readMCIFile.exe: readMCIFile.o
	$(MSVC)/link readMCIFile.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"readMCIFile.exe"

readMCIFile.o: readMCIFile.cpp ../../CUDAMCMLio.h ../../CUDAMCMLio.c
	$(MSVC)/cl readMCIFile.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/I.. \
		/c /Fo"readMCIFile.o"


#####################################################################

detectBoundaryCollision.exe: detectBoundaryCollision.o
	$(MSVC)/link detectBoundaryCollision.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"detectBoundaryCollisions.exe"

detectBoundaryCollision.o: detectBoundaryCollision.cpp geometrylib.c.cpp
	$(MSVC)/cl detectBoundaryCollision.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/I.. \
		/c /Fo"detectBoundaryCollision.o"
		
#TODO geometrylib.c.cpp must be created with cl2cpp.exe from geometrylib.c

clean:
	del *.exe *.o *.preprocessed.c* *.bin.* *.pdb *.ilk
