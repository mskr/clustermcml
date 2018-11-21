// C
//TODO platform-switch to include either direct.h (+redefines) or unistd.h
#include <direct.h> // _getcwd, _stat
#include <string.h> // memcpy
#include <assert.h> // assert

// C++
//TODO get rid of everything that uses exceptions
#include <fstream> // std::ifstream, std::ofstream

// Own
#include "Log.h"
#include "clerr2str.h"
#define DEBUG
#include "clcheck.h" // CL macro
#include "mpicheck.h" // MPI macro

// Interfaces for external code
const char* getCLKernelName();
#ifdef GL_VISUALIZATION
void createGLContexts(void* outDeviceContext, void* outRenderContext);
#endif

// Private module variables
static int rank_ = 0; // MPI process id
static char* src_; // CL kernel source
static cl_device_id* devices_;
static cl_program program_;
static bool initialized_ = false;

/**
*
*/
int getRank() {
	assert(initialized_);
	return rank_;
}

/**
*
*/
void readCLSourceCode(char* kernelfile, char** outsrc, size_t* outlen) {
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
		out << "Could not open kernel file \"" << kernelfile << "\"." << '\n';
		MPI(Abort, MPI_COMM_WORLD, 1);
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

/**
*
*/
void writeCLByteCode(char* kernelfile, cl_program program, cl_uint devicecount, char* devicenames) {
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

/**
*
*/
void broadcastCLSourceCode(char** src, size_t* len) {
	MPI(Bcast, len, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if (rank_ != 0) {
		char* s = new char[*len];
		*src = s;
	}
	MPI(Bcast, *src, (int)*len, MPI_CHAR, 0, MPI_COMM_WORLD);
}

/**
*
*/
void createCLContext(cl_context* outcontext, 
cl_uint* outdevicecount, char* outdevicenames, cl_uint* outplatformcount, char* outplatformnames, size_t memc) {
	// query platform and device info
	cl_platform_id platarr[2]; // detect max 2 platforms
	cl_device_id devarr[5]; // detect max 5 devices
	cl_uint num_platforms_found, num_devices_found;
	CL(GetPlatformIDs, 2, platarr, &num_platforms_found);
	if (outplatformcount) {
		*outplatformcount = num_platforms_found;
	}
	if (outplatformnames && memc > 128 * 2 + 2) {
		size_t platform_name_len = 0, platform_version_len = 0;
		CL(GetPlatformInfo, platarr[0], CL_PLATFORM_NAME, 64, outplatformnames, &platform_name_len);
		memcpy(outplatformnames + platform_name_len - 1, "/", 1); // append cl version after backslash
		CL(GetPlatformInfo, platarr[0], CL_PLATFORM_VERSION, 63, outplatformnames + platform_name_len, &platform_version_len);
		platform_name_len += platform_version_len;
		if (num_platforms_found > 1) {
			memcpy(outplatformnames + platform_name_len - 1, ", ", 2);
			CL(GetPlatformInfo, platarr[1], CL_PLATFORM_NAME, 128, outplatformnames + platform_name_len + 1, NULL);
		}
	}
	// device type: on Radeon5850 clCreateCommandQueue crashed when using all, i.e. shared context for GPU and CPU
	CL(GetDeviceIDs, platarr[0], CL_DEVICE_TYPE_GPU, 5, devarr, &num_devices_found);
	if (outdevicecount) {
		*outdevicecount = num_devices_found;
	}
	if (outdevicenames && memc > 128 * 5 + 8) {
		size_t off = 0;
		for (int i = 0; i < num_devices_found; i++) {
			size_t device_name_len;
			CL(GetDeviceInfo, devarr[i], CL_DEVICE_NAME, 128, outdevicenames + off, &device_name_len);
			off += device_name_len;
			if (i < num_devices_found - 1) {
				memcpy(outdevicenames + (off++) - 1, ", ", 2);
			}
		}
	}
	// create context using all devices of first platform
	#ifdef GL_VISUALIZATION // note that gl makes only sense for single gpu
	cl_context_properties gdiContext, glContext;
	createGLContexts(&gdiContext, &glContext);
	#endif
	cl_context_properties context_platform[] = {
		#ifdef GL_VISUALIZATION
		CL_WGL_HDC_KHR, gdiContext,
		CL_GL_CONTEXT_KHR, glContext,
		#endif
		CL_CONTEXT_PLATFORM, (cl_context_properties)platarr[0], 
		0
	};
	devices_ = new cl_device_id[num_devices_found];
	for (int i = 0; i < num_devices_found; i++) {
		devices_[i] = devarr[i];
	}
	cl_context context = CLCREATE(Context, context_platform, num_devices_found, devices_, NULL, NULL);
	if (outcontext) {
		*outcontext = context;
	}
}

/**
*
*/
void compileCLSourceCode(char* cl_compiler_options, char* src, size_t len,
cl_context context, cl_device_id* devices, cl_uint devicecount, cl_program* outprogram, char* outerr, size_t memc) {
	// To share c header files with kernel append "-I <cwd>" to compiler options
	// (some implementations compile the source as a file in some temp folder)
	int i = 0; for (; cl_compiler_options[i] != '\0'; i++);
	char* opts = (char*)malloc(i + 4 + 256);
	memcpy(opts, cl_compiler_options, i);
	memcpy(&opts[i], " -I ", 4);
	_getcwd(&opts[i + 4], 256);
	const char* programstr[1] = {src};
	size_t programlen[1] = {len};
	cl_program program = CLCREATE(ProgramWithSource, context, 1, programstr, programlen);
	CL(BuildProgram, program, devicecount, devices, opts, NULL, NULL);
	free(opts);
	// Query errors during compilation
	for (int j = 0; j < devicecount; j++) {
		cl_build_status status;
		CL(GetProgramBuildInfo, program, devices[j], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
		if (status == CL_BUILD_ERROR) {
			CL(GetProgramBuildInfo, program, devices[j], CL_PROGRAM_BUILD_LOG, memc, outerr, NULL);
			outerr[0] = 'E';
			break;
		}
	}
	*outprogram = program;
}

/**
*
*/
void createCLCommandQueue(cl_context context, cl_device_id device, cl_command_queue* outcmdqueue) {
	cl_command_queue_properties queueProps = CL_QUEUE_PROFILING_ENABLE/* | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE*/;
	*outcmdqueue = CLCREATE(CommandQueue/*WithProperties*/, context, device, queueProps);
}

/**
*
*/
void usage()  {
	out << "Usage:" << '\n';
	out << "Argument 1: OpenCL kernel file" << '\n';
	out << "Argument 2 (optional): OpenCL compiler options" << '\n';
	out << "Argument 3 (optional): Options for " << getCLKernelName() << '\n';
	MPI(Abort, MPI_COMM_WORLD, 1);
}

/**
*
*/
void initCluster(int nargs, char** args, int* outProcessCount, 
cl_context* outContext, cl_command_queue* outCommandQueue, cl_kernel* outKernel,
size_t* outTotalThreadCount, size_t* outSimdThreadCount) {
	if(nargs < 2) usage();
	char* kernelfile = args[1];
	char* cl_compiler_options = nargs >= 3 ? args[2] : "";
	MPI(Init, &nargs, &args);
	MPI(Comm_rank, MPI_COMM_WORLD, &rank_);

	cl_uint deviceCount, platformCount;
	char deviceNames[1280];
	char platformNames[1280];
	cl_context context;
	createCLContext(&context, &deviceCount, deviceNames, &platformCount, platformNames, 1280);
	if(nargs < 2) usage();
	size_t len;
	if (rank_ == 0) {
		// Only main process should do I/O
		out << platformCount << " CL platforms: " << platformNames << '\n';
		out << deviceCount << " CL devices: " << deviceNames << '\n';
		readCLSourceCode(kernelfile, &src_, &len);
	}
	broadcastCLSourceCode(&src_, &len);
	char compilationErrors[16384];
	compilationErrors[0] = '\0';
	compileCLSourceCode(cl_compiler_options, src_, len, context, devices_, deviceCount, &program_, compilationErrors, 16384);
	if (rank_ == 0) {
		if (compilationErrors[0] != '\0') {
			out << "CL Errors (if you don't see any increase array size and recompile):\n";
			out << compilationErrors << '\n' << Log::flush;
			MPI(Abort, MPI_COMM_WORLD, 1);
		}
		writeCLByteCode(kernelfile, program_, deviceCount, deviceNames);
	}
	createCLCommandQueue(context, devices_[0], outCommandQueue);
	MPI(Comm_size, MPI_COMM_WORLD, outProcessCount);
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
	cl_kernel kernel = CLCREATE(Kernel, program_, getCLKernelName());
	size_t warpSize;
	CL(GetKernelWorkGroupInfo, kernel, devices_[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &warpSize, NULL);
	size_t simdThreadCount;
	CL(GetKernelWorkGroupInfo, kernel, devices_[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &simdThreadCount, NULL);
	simdThreadCount -= simdThreadCount % warpSize;
	cl_uint multiprocessorCount;
	CL(GetDeviceInfo, devices_[0], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &multiprocessorCount, NULL);
	size_t totalThreadCount = multiprocessorCount * simdThreadCount;
	*outContext = context;
	*outKernel = kernel;
	*outSimdThreadCount = simdThreadCount;
	*outTotalThreadCount = totalThreadCount;
	initialized_ = true;
}

/**
*
*/
void cleanupCluster(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel) {
	assert(initialized_);
	CL(ReleaseKernel, kernel);
	CL(ReleaseProgram, program_);
	CL(ReleaseCommandQueue, cmdQueue);
	CL(ReleaseContext, context);
	delete[] src_; // allocated in readCLSourceCode
	delete[] devices_; // allocated in createCLContext
	MPI(Finalize);
	initialized_ = false;
}