/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in 
 *	Think C. If 1, remember to turn on "Generate profiler 
 *	calls" in the options menu. 
 ****/
#define THINKCPROFILER 0	

/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "mcml.h"

/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *);
void FreeData(InputStruct, OutStruct *);
double Rspecular(LayerStruct * );
void LaunchPhoton(double, LayerStruct *, PhotonStruct *);
void HopDropSpin(InputStruct  *,PhotonStruct *,OutStruct *);
void SumScaleResult(InputStruct, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *);


/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on 
 *	screen, return the real time since F=0. 
 *
 *	If F = 2, same as F=1 except no printing.  
 *
 *	Note that clock() and time() return user time and real 
 *	time respectively.
 *	User time is whatever the system allocates to the 
 *	running of the program; 
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *	
 *	clock() only hold 16 bit integer, which is about 32768 
 *	clock ticks.
 ****/
time_t PunchTime(char F, char *Msg)
{
#if GNUCC
  return(0);
#else
  static clock_t ut0;	/* user time reference. */
  static time_t  rt0;	/* real time reference. */
  double secs;
  char s[STRLEN];
  
  if(F==0) {
    ut0 = clock();
    rt0 = time(NULL);
    return(0);
  }
  else if(F==1)  {
    secs = (clock() - ut0)/(double)CLOCKS_PER_SEC;
    if (secs<0) secs=0;	/* clock() can overflow. */
    sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n", 
	    secs, secs/3600.0, Msg);
    puts(s);
    strcpy(Msg, s);
    return(difftime(time(NULL), rt0));
  }
  else if(F==2) return(difftime(time(NULL), rt0));
  else return(0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)	
{
  time_t now, done_time;
  struct tm *date;
  char s[80];
  
  now = time(NULL);
  date = localtime(&now);
  strftime(s, 80, "%H:%M %x", date);
  //printf("Now %s, ", s);
  
  done_time = now + 
			(time_t) (PunchTime(2,"")/(double)P1*(Pt-P1));
  date = localtime(&done_time);
  strftime(s, 80, "%H:%M %x", date);
  //printf("End %s\n", s);
}

/***********************************************************
 *	Report time and write results. 
 ****/
void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
  char time_report[STRLEN];
  
  strcpy(time_report, " Simulation time of this run.");
  PunchTime(1, time_report);

  SumScaleResult(In_Parm, &Out_Parm);
  WriteResult(In_Parm, Out_Parm, time_report);
}

/***********************************************************
 *	Get the file name of the input data file from the 
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc,
					  char * argv[],
					  char * input_filename)
{
	//printf("%s\n", "yolo");
  if(argc>=2) {			/* filename in command line */
	sprintf(input_filename, "%s", argv[1]);
    // strcpy(input_filename, argv[1]);
  }
  else
    input_filename[0] = '\0';
//printf("%s\n", "swag");
} 



    
/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(short NumRuns, InputStruct *In_Ptr)
{
  register long i_photon;	
	/* index to photon. register for speed.*/
  OutStruct out_parm;		/* distribution of photons.*/
  PhotonStruct photon;
  long num_photons = In_Ptr->num_photons, photon_rep=10;

#if THINKCPROFILER
  InitProfile(200,200); cecho2file("prof.rpt",0, stdout);
#endif
    
  InitOutputData(*In_Ptr, &out_parm);
  out_parm.Rsp = Rspecular(In_Ptr->layerspecs);	
  i_photon = num_photons;
  PunchTime(0, "");
    
  do {
    LaunchPhoton(out_parm.Rsp, In_Ptr->layerspecs, &photon);
    int bounces = 0;
    do {
      In_Ptr->bounceCount = bounces;
      HopDropSpin(In_Ptr, &photon, &out_parm);
      bounces++;
    } while (!photon.dead);
    printf("\r%d/%d", num_photons - i_photon + 1, num_photons);
  } while(--i_photon);
  printf("\n");
    
#if THINKCPROFILER
  exit(0);
#endif
    
  ReportResult(*In_Ptr, out_parm);
  FreeData(*In_Ptr, &out_parm);
}

/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/
char main_original_mcml(int argc, char *argv[]) 
{
	//printf("%d\n", STRLEN);
  char input_filename[STRLEN];
  FILE *input_file_ptr;
  short num_runs;	/* number of independent runs. */
  InputStruct in_parm;

  ShowVersion("Version 1.2, 1993");
  GetFnameFromArgv(argc, argv, input_filename);
  input_file_ptr = GetFile(input_filename);
  //printf("%s\n", "hi");
  CheckParm(input_file_ptr, &in_parm);	
  //printf("%s\n", "2");
  num_runs = ReadNumRuns(input_file_ptr);
  //printf("%s\n", "3");
  
  while(num_runs--)  {
    ReadParm(input_file_ptr, &in_parm);
	DoOneRun(num_runs, &in_parm);
  }
  
  fclose(input_file_ptr);
  return(0);
}


#include <assert.h>

#define cl_int int
#define cl_context int
#define cl_command_queue int
#define cl_kernel int
#define cl_mem int
#define cl_event int
#define cl_ulong unsigned long int

void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
  int* inputBufferCount, size_t* inputBufferSizes,
  int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount);
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
  size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);

int main(int nargs, char* args[]) {
  assert(nargs == 2);
  char* filename = args[1];
  {
    int inputBufferCount = 0, outputBufferCount = 0;
    size_t inputBufferSizes[10], outputBufferSizes[10];
    int inputBuffers[10], outputBuffers[10];
    allocCLKernelResources(1, "", filename, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10);
  }

  runCLKernel(0, 0, 0, 0, 0, 1, 1, 1, 0);

  return main_original_mcml(nargs, args);

}