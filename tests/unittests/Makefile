MSVC                = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\bin\Hostx86\x86"
MSVC_INCLUDE        = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\include" # STL
MSVC_INCLUDE_UCRT   = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\ucrt" # crtdefs.h
MSVC_INCLUDE_UM     = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\um"
MSVC_INCLUDE_SHARED = "C:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\shared"
MSVC_LIB            = "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.14.26428\lib\x86"
MSVC_LIB_UCRT       = "C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\ucrt\x86"
MSVC_LIB_UM         = "C:\Program Files (x86)\Windows Kits\8.1\Lib\winv6.3\um\x86"


README.txt:
	echo Available make targets: checkBoundaries.exe, setupInputArrays.exe, readMCIFile.exe, detectBoundaryCollision.exe > README.txt


#####################################################################


checkBoundaries.exe: checkBoundaries.o clmem-windows.o
	$(MSVC)/link checkBoundaries.o clmem-windows.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"checkBoundaries.exe"

checkBoundaries.o: checkBoundaries.cpp ../../CUDAMCMLio.h ../../CUDAMCMLio.c ../../clcheck.h ../../clmem.h ../../geometrylib.h ../../Layer.h ../../Boundary.h
	$(MSVC)/cl checkBoundaries.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/I.. \
		/D"NO_GPU" \
		/c /Fo"checkBoundaries.o"


#####################################################################

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

# Link
detectBoundaryCollision.exe: detectBoundaryCollision.o
	$(MSVC)/link detectBoundaryCollision.o singlemcml-geometrylib-windows.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"detectBoundaryCollision.exe"

# Compile
detectBoundaryCollision.o: detectBoundaryCollision.cpp.transpiled.cpp.preprocessed.cpp singlemcml-geometrylib-windows.o ../../clcheck.h ../../clmem.h ../../geometrylib.h ../../Layer.h ../../Boundary.h
	$(MSVC)/cl detectBoundaryCollision.cpp.transpiled.cpp.preprocessed.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/D"NO_GPU" \
		/c /Fo"detectBoundaryCollision.o"

# Preprocess transpiler stuff (GLM etc.)
detectBoundaryCollision.cpp.transpiled.cpp.preprocessed.cpp: detectBoundaryCollision.cpp.transpiled.cpp
	$(MSVC)/cl detectBoundaryCollision.cpp.transpiled.cpp /c /P \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/D"NO_GPU" \
		/P /Fi"detectBoundaryCollision.cpp.transpiled.cpp.preprocessed.cpp"


# Transpile
detectBoundaryCollision.cpp.transpiled.cpp: detectBoundaryCollision.cpp cl2cpp.exe
	cl2cpp detectBoundaryCollision.cpp -v

# Preprocess
# detectBoundaryCollision.preprocessed.cpp: detectBoundaryCollision.cpp
# 	$(MSVC)/cl detectBoundaryCollision.cpp /c /P \
# 		/I$(MSVC_INCLUDE) \
# 		/I$(MSVC_INCLUDE_UCRT) \
# 		/I"../.." \
# 		/D"NO_GPU" \
# 		/P /Fi"detectBoundaryCollision.preprocessed.cpp"
		

# 5. Now we can compile the kernel as C++
singlemcml-geometrylib-windows.o: singlemcml-geometrylib-windows.preprocessed.cpp
	$(MSVC)/cl singlemcml-geometrylib-windows.preprocessed.cpp /c /Zi /W3 /Fo"singlemcml-geometrylib-windows.o"

# 4. Another preprocessing step
singlemcml-geometrylib-windows.preprocessed.cpp: geometrylib.preprocessed.c.transpiled.cpp
	$(MSVC)/cl geometrylib.preprocessed.c.transpiled.cpp /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/I"../.." \
		/P /Fi"singlemcml-geometrylib-windows.preprocessed.cpp"

# 3. Transpile the preprocessed cl code to cpp
geometrylib.preprocessed.c.transpiled.cpp: geometrylib.preprocessed.c cl2cpp.exe
	cl2cpp geometrylib.preprocessed.c -v

# 2. Preprocess file to resolve includes
geometrylib.preprocessed.c: ../../geometrylib.c ../../geometrylib.h
	$(MSVC)/cl "../../geometrylib.c" /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/D"NO_GPU" \
		/P /Fi"geometrylib.preprocessed.c"

# 1. Make transpiler if not present
cl2cpp.exe: cl2cpp.o
	$(MSVC)/link cl2cpp.o /DEBUG \
		/LIBPATH:$(MSVC_LIB) \
		/LIBPATH:$(MSVC_LIB_UCRT) \
		/LIBPATH:$(MSVC_LIB_UM) \
		/OUT:"cl2cpp.exe"

cl2cpp.o: ../../cl2cpp.cpp
	$(MSVC)/cl "../../cl2cpp.cpp" /c /Zi /W3 \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/c /Fo"cl2cpp.o"


#####################################################################


clmem-windows.o: clmem.preprocessed.cpp
	$(MSVC)/cl clmem.preprocessed.cpp /c /W3 /Fo"clmem-windows.o"

clmem.preprocessed.cpp: ../../clmem.cpp
	$(MSVC)/cl "../../clmem.cpp" /c \
		/I$(MSVC_INCLUDE) \
		/I$(MSVC_INCLUDE_UCRT) \
		/D"NO_GPU" \
		/P /Fi"clmem.preprocessed.cpp"

clean:
	del *.exe *.o *.preprocessed.c* *.bin.* *.pdb *.ilk *.transpiled.cpp