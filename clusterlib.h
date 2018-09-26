/**
*
*/
int getRank();

/**
*
*/
void readCLSourceCode(char* kernelfile, char** outsrc, size_t* outlen);

/**
*
*/
void writeCLByteCode(char* kernelfile, cl_program program, cl_uint devicecount, char* devicenames);

/**
*
*/
void broadcastCLSourceCode(char** src, size_t* len);

/**
*
*/
void createCLContext(cl_context* outcontext, cl_device_id** outdevices, 
cl_uint* outdevicecount, char* outdevicenames, cl_uint* outplatformcount, char* outplatformnames, size_t memc);

/**
*
*/
void compileCLSourceCode(char* cl_compiler_options, char* src, size_t len,
cl_context context, cl_device_id* devices, cl_uint devicecount, cl_program* outprogram, char* outerr, size_t memc);

/**
*
*/
void createCLCommandQueue(cl_context context, cl_device_id device, cl_command_queue* outcmdqueue);

/**
*
*/
void initCluster(int nargs, char** args, int* outProcessCount, 
cl_context* outContext, cl_command_queue* outCommandQueue, cl_kernel* outKernel,
size_t* outTotalThreadCount, size_t* outSimdThreadCount);

/**
*
*/
void cleanupCluster(int inputBufferCount, cl_mem* inputBuffers, int outputBufferCount, cl_mem* outputBuffers,
cl_context context, cl_command_queue cmdQueue, cl_kernel kernel);


//TODO manage better work distribution:
struct CLDeviceInfo {
	cl_command_queue queue;
	int totalThreadCount;
	int simdThreadCount;
};
struct CLContextInfo {
	int deviceCount;
	CLDeviceInfo devices[10];
};
struct ClusterInfo {
	int processCount;
	int clContextCount;
	CLContextInfo clContexts[10];
};