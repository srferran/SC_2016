#include <stdio.h>
#include <iostream>
#include "CImg.h" 
#include <CL/cl.h>




// Source will be read from disk
#define MAX_SOURCE_SIZE (0x100000)

using namespace std;
using namespace cimg_library;

CImg<float> pattern_matching(CImg<unsigned char> &img, CImg<unsigned char> &pat) 
{
  // Output. We assume input image has a border of 15 pixels
  CImg<float> out(img.width()-15,img.height()-15,1,1,0);

  // Kernel source
  FILE *fp;
  char *programSource;
  size_t programSize;

  // For the timer 
  struct timeval start, end;
  float mtime, seconds, useconds;


  fp = fopen("kernel.cl", "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  programSource = (char *) malloc(sizeof(char) * MAX_SOURCE_SIZE);
  programSize = fread(programSource, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  cl_int status;  

  //-----------------------------------------------------
  // Discover and initialize the platforms
  //-----------------------------------------------------

  cl_uint numPlatforms = 0;
  cl_platform_id *platforms = NULL;

  status = clGetPlatformIDs(0, NULL, &numPlatforms);

  platforms =   
    (cl_platform_id*)malloc(
	numPlatforms*sizeof(cl_platform_id));

  status = clGetPlatformIDs(numPlatforms, platforms, 
      NULL);

  //-----------------------------------------------------
  // Discover and initialize the devices
  //----------------------------------------------------- 

  cl_uint numDevices = 0;
  cl_device_id *devices = NULL;

  status = clGetDeviceIDs(
      platforms[0], 
      CL_DEVICE_TYPE_ALL, 
      0, 
      NULL, 
      &numDevices);

  devices = 
    (cl_device_id*)malloc(
	numDevices*sizeof(cl_device_id));

  status = clGetDeviceIDs(
      platforms[0], 
      CL_DEVICE_TYPE_ALL,        
      numDevices, 
      devices, 
      NULL);

  //-----------------------------------------------------
  // Create a context
  //----------------------------------------------------- 

  cl_context context = NULL;

  context = clCreateContext(
      NULL, 
      numDevices, 
      devices, 
      NULL, 
      NULL, 
      &status);

  if (status != CL_SUCCESS) {
    printf("ERROR: creating context.\n");
    exit(1);
  }

  //-----------------------------------------------------
  // Create a command queue
  //----------------------------------------------------- 

      cl_command_queue cmdQueue;

      cmdQueue = clCreateCommandQueue(
      context, 
      devices[0], 
      CL_QUEUE_PROFILING_ENABLE, 
      &status);
  
  //-----------------------------------------------------
  // Create and compile the program
  //----------------------------------------------------- 

  cl_program program = clCreateProgramWithSource(
      context, 
      1, 
      (const char **) &programSource,                                 
      (const size_t *) &programSize, 
      &status);

  if (status != CL_SUCCESS) {
    printf("ERROR: creating program.\n");
    exit(1);
  }

  status = clBuildProgram(
      program, 
      numDevices, 
      devices, 
      NULL, 
      NULL, 
      NULL);

  if (status != CL_SUCCESS) {
    cl_build_status build_status;
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);

    char *build_log;
    size_t ret_val_size;
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);

    build_log = (char *) malloc(ret_val_size+1);
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
    build_log[ret_val_size] = '\0';
    printf("BUILD LOG: \n %s", build_log);
    free(build_log);
    exit(1);
  }

  //-----------------------------------------------------
  // We start the timer 
  //----------------------------------------------------- 

  gettimeofday(&start, NULL);

  //-----------------------------------------------------
  // Create device buffers
  //----------------------------------------------------- 

  cl_mem bufferImg; 
  cl_mem bufferPat; 
  cl_mem bufferOut; 

  bufferImg = clCreateBuffer(
      context, 
      CL_MEM_READ_ONLY,                         
      sizeof(char) * img.height() * img.width(), 
      NULL, 
      &status);

  bufferPat = clCreateBuffer(
      context, 
      CL_MEM_READ_ONLY,                         
      sizeof(char) * pat.height() * pat.width(), 
      NULL, 
      &status);

  bufferOut = clCreateBuffer(
      context, 
      CL_MEM_WRITE_ONLY,                 
      sizeof(float) * out.height() * out.width(), 
      NULL, 
      &status);

  //-----------------------------------------------------
  // Write host data to device buffers
  //----------------------------------------------------- 

  status = clEnqueueWriteBuffer(
      cmdQueue, 
      bufferImg, 
      CL_FALSE, 
      0, 
      sizeof(char) * img.height() * img.width(),                         
      img._data, 
      0, 
      NULL, 
      NULL);

  status = clEnqueueWriteBuffer(
      cmdQueue, 
      bufferPat, 
      CL_FALSE, 
      0, 
      sizeof(char) * pat.height() * pat.width(),                                  
      pat._data, 
      0, 
      NULL, 
      NULL);

  //-----------------------------------------------------
  // Create the kernel
  //----------------------------------------------------- 

  cl_kernel kernel = NULL;

  kernel = clCreateKernel(program, "pattern_matching", &status);

  if (status != CL_SUCCESS) {
    printf("ERROR: creating kernel.\n");
    exit(1);
  }

  //-----------------------------------------------------
  // Set the kernel arguments
  //----------------------------------------------------- 

  int height = img.height();
  int width  = img.width();

  status  = clSetKernelArg(
      kernel, 
      0, 
      sizeof(cl_mem), 
      &bufferImg);

  status |= clSetKernelArg(
      kernel, 
      1, 
      sizeof(cl_mem), 
      &bufferPat);

  status |= clSetKernelArg(
      kernel, 
      2, 
      sizeof(int), 
      &height);

  status |= clSetKernelArg(
      kernel, 
      3, 
      sizeof(int), 
      &width);

  status |= clSetKernelArg(
      kernel, 
      4, 
      sizeof(cl_mem), 
      &bufferOut);

  if (status != CL_SUCCESS) {
    printf("ERROR: setting arguments.\n");
    exit(1);
  }

  clFinish(cmdQueue);

  //-----------------------------------------------------
  // Configure the work-item structure
  //----------------------------------------------------- 
  
  size_t localWorkSize[] = {32, 32};
  size_t globalWorkSize[] = {256,256};
  
  //-----------------------------------------------------
  // Enqueue the kernel for execution
  //----------------------------------------------------- 

  cl_event event;
  cl_ulong time_start, time_end;
  double total_time;

  status = clEnqueueNDRangeKernel(
      cmdQueue, 
      kernel, 
      2, 
      NULL, 
      globalWorkSize, 
      localWorkSize, 
      0, 
      NULL, 
      &event);

  clFinish(cmdQueue);

  clGetEventProfilingInfo(event,  
      CL_PROFILING_COMMAND_START, 
      sizeof(time_start), 
      &time_start, 
      NULL);

  clGetEventProfilingInfo(event, 
      CL_PROFILING_COMMAND_END, 
      sizeof(time_end), 
      &time_end, 
      NULL);

  total_time = time_end - time_start;
  printf("Execution of kernel time at GPU = %0.5f ms\n", total_time / 1000000.0);
  printf("\n");
  //-----------------------------------------------------
  // Read the output buffer back to the host
  //----------------------------------------------------- 

  clEnqueueReadBuffer(
      cmdQueue, 
      bufferOut, 
      CL_TRUE, 
      0, 
      sizeof(float) * out.height() * out.width(), 
      out._data, 
      0, 
      NULL, 
      NULL);

  //-----------------------------------------------------
  // Get again timer 
  //----------------------------------------------------- 

  gettimeofday(&end, NULL);	
  seconds = end.tv_sec - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;

  mtime = ((seconds) * 1000.0 + useconds / 1000.0);
  //mostrem el temps emprat
  printf("Elapsed time of pattmatch function = %0.5f ms", mtime);

  //-----------------------------------------------------
  // Release OpenCL resources
  //----------------------------------------------------- 

  // Free OpenCL resources
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(cmdQueue);
  clReleaseMemObject(bufferImg);
  clReleaseMemObject(bufferPat);
  clReleaseMemObject(bufferOut);
  clReleaseContext(context);

  free(programSource);
  free(platforms);
  free(devices);

  //--------------------------------------------------
  // Return ouput
  //--------------------------------------------------

  return out;
}
