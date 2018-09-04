// C
//TODO platform-switch to include either direct.h (+redefines) or unistd.h
#include <direct.h> // _getcwd, _stat
#include <string.h> // memcpy

// C++
#include <iostream> // std::cout
#include <fstream> // std::ifstream, std::ofstream

// Own
#include "clerr2str.c"
#define DEBUG
#include "clcheck.h"

// MPI process id
static int rank = 0;

// interfaces for external code
const char* getCLKernelName();
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
	int* inputBufferCount, size_t* inputBufferSizes,
	int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount);
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
	size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);
void createGLContexts(void* outDeviceContext, void* outRenderContext);
void createGLBuffer(size_t size, void* outBuffer);
void runGLRenderLoop();

void usage()  {
	std::cout << "Usage:" << std::endl;
	std::cout << "Argument 1: OpenCL kernel file" << std::endl;
	std::cout << "Argument 2 (optional): OpenCL compiler options" << std::endl;
	std::cout << "Argument 3 (optional): Options for " << getCLKernelName() << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
}

void readCLSourceCode(int nargs, char** args, char** outsrc, size_t* outlen) {
	if(nargs < 2) usage();
	char kfilepath[256];
	_getcwd(kfilepath, 256);
	int stri = 0;
	char c = kfilepath[0];
	while (c != '\0') {
		stri++;
		c = kfilepath[stri];
	}
	kfilepath[stri] = '/';
	stri++;
	char* kernelfile = args[1];
	c = kernelfile[0];
	int j = 0;
	while (c != '\0') {
		kfilepath[stri] = c;
		stri++;
		j++;
		c = kernelfile[j];
	}
	kfilepath[stri] = '\0';
	std::ifstream kfilestream(kfilepath, std::fstream::in | std::ifstream::binary);
	if (!kfilestream.is_open()) {
		std::cout << "Could not open kernel file \"" << kernelfile << "\"." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	size_t kernellen = 0;
	struct _stat kfilestat;
	_stat(kfilepath, &kfilestat);
	kernellen = kfilestat.st_size;
	char* kernelstr = new char[kernellen + 1]();
	kfilestream.read(kernelstr, kernellen);
	kernelstr[kernellen] = '\0';
	*outsrc = kernelstr;
	*outlen = kernellen + 1;
}

void writeCLByteCode(int nargs, char** args, cl_program program, cl_uint devicecount, char* devicenames) {
	if(nargs < 2) usage();
	char kfilepath[256];
	_getcwd(kfilepath, 256);
	int stri = 0;
	char c = kfilepath[0];
	while (c != '\0') {
		stri++;
		c = kfilepath[stri];
	}
	kfilepath[stri] = '/';
	stri++;
	char* kernelfile = args[1];
	c = kernelfile[0];
	int j = 0;
	while (c != '\0') {
		kfilepath[stri] = c;
		stri++;
		j++;
		c = kernelfile[j];
	}
	kfilepath[stri] = '\0';
	size_t* binsizes = new size_t[devicecount];
	CL(GetProgramInfo, program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t) * devicecount, binsizes, NULL);
	char** kernelbins = new char*[devicecount];
	size_t totalbinsizes = 0;
	for (int i = 0; i < devicecount; i++) {
		if (binsizes[i] > 0) {
			kernelbins[i] = new char[binsizes[i]]();
			totalbinsizes += binsizes[i];
		}
	}
	if (totalbinsizes == 0) {
		delete[] binsizes;
		delete[] kernelbins;
		return; // no binaries available
	}
	CL(GetProgramInfo, program, CL_PROGRAM_BINARIES, devicecount * sizeof(char*), kernelbins, NULL);
	char** device_names = new char*[devicecount];
	int idx = 0;
	int len = 0;
	int device_idx = 0;
	int* device_name_lengths = new int[devicecount];
	while (devicenames[idx] != '\0') {
		if (devicenames[idx] == ',' && devicenames[idx + 1] == ' ') {
			device_names[device_idx] = new char[len + 1];
			memcpy(device_names[device_idx], &devicenames[idx - len], len);
			device_names[device_idx][len] = '\0';
			device_name_lengths[device_idx] = len + 1;
			len = 0;
			device_idx++;
		} else if (!(devicenames[idx - 1] == ',' && devicenames[idx] == ' ')) {
			len++;
		}
		idx++;
	}
	device_names[device_idx] = new char[len + 1];
	memcpy(device_names[device_idx], &devicenames[idx - len], len);
	device_names[device_idx][len] = '\0';
	device_name_lengths[device_idx] = len + 1;
	for (int i = 0; i < devicecount; i++) {
		memcpy(&kfilepath[stri], ".bin.", 5);
		memcpy(&kfilepath[stri + 5], device_names[i], device_name_lengths[i]);
		std::ofstream outfile(kfilepath);
		outfile.write(kernelbins[i], binsizes[i]);
		delete[] kernelbins[i];
		delete[] device_names[i];
	}
	delete[] device_name_lengths;
	delete[] device_names;
	delete[] binsizes;
	delete[] kernelbins;
}

void broadcastCLSourceCode(char** src, size_t* len) {
	MPI_Bcast(len, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if (rank != 0) {
		char* s = new char[*len];
		*src = s;
	}
	MPI_Bcast(*src, (int)*len, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void createCLContext(cl_context* outcontext, cl_device_id** outdevices, 
cl_uint* outdevicecount, char* outdevicenames, cl_uint* outplatformcount, char* outplatformnames, size_t memc) {
	// query platform and device info
	cl_platform_id platforms[2]; // detect max 2 platforms
	cl_device_id devices[5]; // detect max 5 devices
	cl_uint num_platforms_found, num_devices_found;
	CL(GetPlatformIDs, 2, platforms, &num_platforms_found);
	if (outplatformcount) {
		*outplatformcount = num_platforms_found;
	}
	if (outplatformnames && memc > 128 * 2 + 2) {
		size_t platform_name_len = 0, platform_version_len = 0;
		CL(GetPlatformInfo, platforms[0], CL_PLATFORM_NAME, 64, outplatformnames, &platform_name_len);
		memcpy(outplatformnames + platform_name_len - 1, "/", 1); // append cl version after backslash
		CL(GetPlatformInfo, platforms[0], CL_PLATFORM_VERSION, 63, outplatformnames + platform_name_len, &platform_version_len);
		platform_name_len += platform_version_len;
		if (num_platforms_found > 1) {
			memcpy(outplatformnames + platform_name_len - 1, ", ", 2);
			CL(GetPlatformInfo, platforms[1], CL_PLATFORM_NAME, 128, outplatformnames + platform_name_len + 1, NULL);
		}
	}
	// device type: on Radeon5850 clCreateCommandQueue crashed when using all, i.e. shared context for GPU and CPU
	CL(GetDeviceIDs, platforms[0], CL_DEVICE_TYPE_GPU, 5, devices, &num_devices_found);
	if (outdevicecount) {
		*outdevicecount = num_devices_found;
	}
	if (outdevicenames && memc > 128 * 5 + 8) {
		size_t off = 0;
		for (int i = 0; i < num_devices_found; i++) {
			size_t device_name_len;
			CL(GetDeviceInfo, devices[i], CL_DEVICE_NAME, 128, outdevicenames + off, &device_name_len);
			off += device_name_len;
			if (i < num_devices_found - 1) {
				memcpy(outdevicenames + (off++) - 1, ", ", 2);
			}
		}
	}
	// create context using all devices of first platform
	#ifdef GL_VISUALIZATION
	cl_context_properties gdiContext, glContext;
	createGLContexts(&gdiContext, &glContext);
	#endif
	cl_context_properties context_platform[] = {
		#ifdef GL_VISUALIZATION
		CL_WGL_HDC_KHR, gdiContext,
		CL_GL_CONTEXT_KHR, glContext,
		#endif
		CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[0], 
		0
	};
	cl_device_id* context_devices = new cl_device_id[num_devices_found];
	for (int i = 0; i < num_devices_found; i++) {
		context_devices[i] = devices[i];
	}
	cl_context context = CLCREATE(Context, context_platform, num_devices_found, context_devices, NULL, NULL);
	if (outdevices) {
		*outdevices = context_devices;
	} else {
		delete[] context_devices;
	}
	if (outcontext) {
		*outcontext = context;
	}
}

void compileCLSourceCode(int nargs, char** args, char* src, size_t len,
cl_context context, cl_device_id* devices, cl_uint devicecount, cl_program* outprogram, char* outerr, size_t memc) {
	char* cl_compiler_options = nargs >= 3 ? args[2] : "";
	//TODO to share c header files with kernel append -I <getcwd> to compiler options
	const char* programstr[1] = {src};
	size_t programlen[1] = {len};
	cl_program program = CLCREATE(ProgramWithSource, context, 1, programstr, programlen);
	CL(BuildProgram, program, devicecount, devices, cl_compiler_options, NULL, NULL);
	// query errors during compilation
	for (int i = 0; i < devicecount; i++) {
		cl_build_status status;
		CL(GetProgramBuildInfo, program, devices[i], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
		if (status == CL_BUILD_ERROR) {
			CL(GetProgramBuildInfo, program, devices[i], CL_PROGRAM_BUILD_LOG, memc, outerr, NULL);
			outerr[0] = 'E';
			break;
		}
	}
	*outprogram = program;
}

void createCLCommandQueue(cl_context context, cl_device_id device, cl_command_queue* outcmdqueue) {
	cl_command_queue_properties queueProps = CL_QUEUE_PROFILING_ENABLE/* | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE*/;
	*outcmdqueue = CLCREATE(CommandQueue/*WithProperties*/, context, device, queueProps);
}

int main(int nargs, char** args) {
	MPI_Init(&nargs, &args);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cl_context context;
	cl_device_id* devices;
	cl_uint deviceCount, platformCount;
	char deviceNames[1280];
	char platformNames[1280];
	createCLContext(&context, &devices, &deviceCount, deviceNames, &platformCount, platformNames, 1280);
	char* src; size_t len;
	if (rank == 0) {
		// Only main process should do I/O
		std::cout << platformCount << " CL platforms: " << platformNames << std::endl;
		std::cout << deviceCount << " CL devices: " << deviceNames << std::endl;
		readCLSourceCode(nargs, args, &src, &len);
	}
	broadcastCLSourceCode(&src, &len);
	cl_program program;
	char compilationErrors[4096];
	compilationErrors[0] = '\0';
	compileCLSourceCode(nargs, args, src, len, context, devices, deviceCount, &program, compilationErrors, 4096);
	if (rank == 0) {
		if (compilationErrors[0] != '\0') {
			std::cout << compilationErrors << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		writeCLByteCode(nargs, args, program, deviceCount, deviceNames);
	}
	cl_command_queue cmdQueue1;
	createCLCommandQueue(context, devices[0], &cmdQueue1);
	int processCount;
	MPI_Comm_size(MPI_COMM_WORLD, &processCount);
	/*****************************************************************************************************************
	// some of the gpu threads share the same register memory and L1 cache, forming "multiprocessors"
	// CL_DEVICE_MAX_COMPUTE_UNITS == number of multiprocessors
	// work groups are definitely inside a multiprocessor
	// the more registers a kernel uses, the smaller a work group should be (less threads)
	// CL_DEVICE_MAX_WORK_GROUP_SIZE == maximum threads in a multiprocessor
	// CL_KERNEL_WORK_GROUP_SIZE == threads of a kernel that fit in a multiprocessor
	// there is no guarantee that the threads will run in parallel
	// Nvidia warps and AMD wavefronts are the subsets of a multiprocessor that share the same instruction dispatcher
	// CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE == warp/wavefront size
	// 
	//                                              | ATI Radeon 5850   | Intel Graphics 4600   | NVIDIA GTX 1080   |
	// --------------------------------------------------------------------------------------------------------------
	// CL_DEVICE_MAX_COMPUTE_UNITS                  | 18                |
	// CL_DEVICE_MAX_WORK_GROUP_SIZE                | 256               |
	// CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE | 64                |
	*****************************************************************************************************************/
	cl_kernel kernel = CLCREATE(Kernel, program, getCLKernelName());
	size_t warpSize;
	CL(GetKernelWorkGroupInfo, kernel, devices[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &warpSize, NULL);
	size_t simdThreadCount;
	CL(GetKernelWorkGroupInfo, kernel, devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &simdThreadCount, NULL);
	simdThreadCount -= simdThreadCount % warpSize;
	cl_uint multiprocessorCount;
	CL(GetDeviceInfo, devices[0], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &multiprocessorCount, NULL);
	size_t totalThreadCount = multiprocessorCount * simdThreadCount;

	char kernelOptions[512];
	CL(GetProgramBuildInfo, program, devices[0], CL_PROGRAM_BUILD_OPTIONS, 512, kernelOptions, NULL);
	char* otherOptions = nargs >= 4 ? args[3] : "";
	int inputBufferCount = 0, outputBufferCount = 0;
	size_t inputBufferSizes[10], outputBufferSizes[10];
	cl_mem inputBuffers[10], outputBuffers[10];
	allocCLKernelResources(totalThreadCount, kernelOptions, otherOptions, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10);
	// Host buffers were created and sizes reported, now create CL buffers
	for (int i = 0; i < inputBufferCount; i++) {
		inputBuffers[i] = CLCREATE(Buffer, context, CL_MEM_READ_ONLY, inputBufferSizes[i], NULL);
	}
	for (int i = 0; i < outputBufferCount; i++) {
		#ifdef GL_VISUALIZATION
		unsigned int glBuf = 0;
		createGLBuffer(outputBufferSizes[i], &glBuf);
		outputBuffers[i] = CLCREATE(FromGLBuffer, context, CL_MEM_READ_WRITE, glBuf);
		#else
		outputBuffers[i] = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, outputBufferSizes[i], NULL);
		#endif
	}
	#ifdef GL_VISUALIZATION
	// glFinish(); no need to finish gl commands as we havent even started the render loop
	for (int i = 0; i < outputBufferCount; i++) {
		CL(EnqueueAcquireGLObjects, cmdQueue1, 1, &outputBuffers[i], 0, NULL, NULL);
	}
	#endif
	runCLKernel(context, cmdQueue1, kernel, inputBuffers, outputBuffers, totalThreadCount, simdThreadCount, processCount, rank);
	#ifdef GL_VISUALIZATION
	for (int i = 0; i < outputBufferCount; i++) {
		CL(EnqueueReleaseGLObjects, cmdQueue1, 1, &outputBuffers[i], 0, NULL, NULL);
	}
	runGLRenderLoop();
	#endif

	for (int i = 0; i < inputBufferCount; i++) {
		CL(ReleaseMemObject, inputBuffers[i]);
	}
	for (int i = 0; i < outputBufferCount; i++) {
		CL(ReleaseMemObject, outputBuffers[i]);
	}
	CL(ReleaseKernel, kernel);
	CL(ReleaseProgram, program);
	CL(ReleaseCommandQueue, cmdQueue1);
	CL(ReleaseContext, context);

	delete[] src;
	delete[] devices;
	MPI_Finalize();

	// hello world example
	// http://phycomp.technion.ac.il/~comphy/classfiles/hellocl.c

	// nvidia optimization: occupancy, ...
	// https://www.cs.cmu.edu/afs/cs/academic/class/15668-s11/www/cuda-doc/OpenCL_Best_Practices_Guide.pdf

	// Single kernel is for data-parallel problems, running multiple kernels in parallel is for task-parallel problems
	// Can use multiple queues or out-of-order-mode to enqueue multiple kernels
	// https://us.fixstars.com/opencl/book/OpenCLProgrammingBook/opencl-programming-practice/

	// work can be further spread on multiple devices using multiple command queues

	// For also using CPU there are options:
	// 1. create a separate cl context and compile kernel for CPU
	// 2. compile kernel with C compiler (macro defines, no vector types) 
	// 3. Automatic scheduling of the submitted commands on OpenCL devices
	//    http://starpu.gforge.inria.fr/doc/html/SOCLOpenclExtensions.html

	// Hot reloading hardcoded parameters
	// http://tuxedolabs.blogspot.com/2018/03/hot-reloading-hardcoded-parameters.html?m=1
	// File watcher: https://facebook.github.io/watchman/
}