MSVC = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128\bin\Hostx86\x86"
MSVC_INCLUDE = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128\include" # STL
MSVC_INCLUDE_UCRT = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\ucrt" # crtdefs.h
MSVC_INCLUDE_UM = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\um"
MSVC_INCLUDE_SHARED = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.16299.0\shared"
MSVC_LIB = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.13.26128/lib/x86"
MSVC_LIB_UCRT = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.16299.0\ucrt\x86"
MSVC_LIB_UM = "C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86"

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
