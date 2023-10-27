#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <sys/time.h>

#include <cmath>
#include <cstring>

#ifdef __APPLE_CC__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

using namespace std;

namespace xppc{
#define LMAX 80      // number of dust loggers
#define LYRS 172     // number of depth points

#define CTX  11
#define CTY  10
#define DIR1 9.3f
#define DIR2 129.3f

#define CX   21
#define CY   19
#define NSTR 240

#define OVER   128   // minimum number of photons per particle (optimized)
#define NPHO   1024  // maximum number of photons propagated by one thread

#define WNUM   512   // number of wavelength slices
#define MAXLYS 172   // maximum number of ice layers
#define MAXGEO 16384 // maximum number of OMs
#define MAXRND 131072 // max. number of random number multipliers

#define XXX 1.e-5f
#define FPI 3.141592653589793f
#define OMR 0.16510f  // DOM radius [m]

#define BUG_with_Unload

#include "pro.cxx"
#include "ini.cxx"

#define SHRT
#include "pro.cxx"
#undef SHRT

  void initialize(float enh = 1.f){ m.set(); q.eff*=enh; }

  unsigned int pmax, pmxo, pn, pk, hquo;

  void setq(){
    char * HQUO=getenv("HQUO");
    hquo=HQUO==NULL?1:atoi(HQUO);
    cerr<<"HQUO(photons/max number of hits)="<<hquo<<endl;
  }

  bool xgpu=false;

  void checkError(cl_int result){
    if(result!=CL_SUCCESS){
      cerr<<"OpenCL Error: "<<result<<endl;
      exit(2);
    }
  }

  vector< pair<cl_platform_id,cl_device_id> > all;

  struct gpu{
    dats d;

    int device, mult;
    cl_uint nblk;
    size_t nthr, ntot;
    unsigned int npho, pmax, pmxo;

    unsigned int old, num;

    float deviceTime;
    cl_event event;

    cl_platform_id pfID;
    cl_device_id devID;
    cl_context ctx;
    cl_command_queue cq;
    cl_program program;
    cl_kernel clkernel;

    cl_mem ed, ez, eh, ep, bf, eo; // pointers to structures on device

    string& replace(string& in, string old, string str){
      string clx;
      int k=0, m=0;
      while((m=in.find(old, k))!=-1){
	clx.append(in, k, m-k);
	k=m+old.length();
	clx.append(str);
      }
      clx.append(in, k, in.length()-k);
      return in=clx;
    }

    gpu(int device) : mult(1), npho(NPHO), old(0), deviceTime(0){
      this->device=device;

      {
	ostringstream o; o<<"NPHO_"<<device;
	char * nph=getenv(o.str().c_str());
	if(nph==NULL) nph=getenv("NPHO");
	if(nph!=NULL) if(*nph!=0){
	    npho=max(0, atoi(nph));
	    cerr<<"Setting NPHO="<<npho<<endl;
	    if(npho<=0){
	      cerr<<"Not using device # "<<device<<"!"<<endl;
	      return;
	    }
	  }
      }

      pfID  = all[device].first; 
      devID = all[device].second;

      string opts;

      cl_device_type dtype;
      clGetDeviceInfo(devID, CL_DEVICE_TYPE, sizeof(cl_device_type), &dtype, NULL);

      {
	cl_int vid;
	clGetDeviceInfo(devID, CL_DEVICE_VENDOR_ID, sizeof(cl_uint), &vid, NULL);
	if(dtype==CL_DEVICE_TYPE_GPU){
	  opts+=" -cl-fast-relaxed-math";
	  switch(vid){
	  case 0x1002: mult=8; break;
	  case 0x10de:
	    opts+=" -cl-nv-maxrregcount=64";
	    opts+=" -cl-nv-verbose"; break;
	  }
	}

	size_t siz;
	clGetDeviceInfo(devID, CL_DEVICE_VENDOR, 0, NULL, &siz);
	char nam[siz];
	clGetDeviceInfo(devID, CL_DEVICE_VENDOR, siz, nam, NULL);
	cerr<<"Device vendor ID=0x"<<hex<<vid<<dec<<" ("<<nam<<") mult="<<mult<<endl;
      }

      {
	ostringstream o; o<<"XMLT_"<<device;
	char * mlt=getenv(o.str().c_str());
	if(mlt==NULL) mlt=getenv("XMLT");
	if(mlt!=NULL) if(*mlt!=0){
	    int aux=atoi(mlt);
	    if(aux>0){
	      mult=aux;
	      cerr<<"Setting XMLT="<<mult<<endl;
	    }
	  }
      }

      {
	cl_int err;

	const cl_context_properties prop[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) pfID, 0};
	ctx = clCreateContext(prop, 1, &devID, NULL, NULL, &err); checkError(err);
	cq = clCreateCommandQueue(ctx, devID, CL_QUEUE_PROFILING_ENABLE, &err); checkError(err);

	string tmp = kernel_source.substr(1, kernel_source.length()-2);
	string source(replace(replace(tmp, "; ", ";\n"), "cl_", ""));

	const char *src = source.c_str();
	size_t length=source.length();

	program = clCreateProgramWithSource(ctx, 1, &src, &length, &err); checkError(err);

	clBuildProgram(program, 1, &devID, opts.c_str(), NULL, NULL);

#ifdef BUG_with_Unload
	// segfaults in some enviroments
#elif CL_VERSION_1_2
	checkError(clUnloadPlatformCompiler(pfID));
#else
	checkError(clUnloadCompiler());
#endif

	cl_build_status status;
	checkError(clGetProgramBuildInfo(program, devID, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL));

	size_t siz;
	checkError(clGetProgramBuildInfo(program, devID, CL_PROGRAM_BUILD_LOG, 0, NULL, &siz));
	char log[siz+1]; log[siz] = '\0';
	checkError(clGetProgramBuildInfo(program, devID, CL_PROGRAM_BUILD_LOG, siz, log, NULL));
	if(strlen(log)>1 || status!=CL_SUCCESS) fprintf(stderr, "BUILD LOG:\n%s\n", log);
	checkError(status);

	clkernel = clCreateKernel(program, "propagate", &err); checkError(err);
      }

      {
	clGetDeviceInfo(devID, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &nblk, NULL);
	clGetKernelWorkGroupInfo(clkernel, devID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &nthr, NULL);

	cl_ulong lmem;
	clGetKernelWorkGroupInfo(clkernel, devID, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &lmem, NULL);

	cerr<<"Running on "<<nblk<<" MPs x "<<nthr<<" threads; Kernel uses: l="<<lmem<<endl;
      }

      if(dtype!=CL_DEVICE_TYPE_GPU) nthr=1;  // appears necessary on CPUs and Intel Phi
      else{
	unsigned int num=0;
	checkError(clSetKernelArg(clkernel, 0, sizeof(unsigned int), &num));
	for(int n=1; n<=6; n++) checkError(clSetKernelArg(clkernel, n, sizeof(cl_mem), NULL));

	do nthr++; while(CL_SUCCESS==clEnqueueNDRangeKernel(cq, clkernel, 1, NULL, &nthr, &nthr, 0, NULL, &event));
	nthr--;

	checkError(clFlush(cq));
      }
      {
	cerr<<"Running on "<<nblk<<" MPs x "<<nthr<<" threads; (corrected/verified)"<<endl;
      }
    }

    void ini(){
      rs_ini();
      d=xppc::d;
      d.gdev=device;

      nblk*=mult;
      ntot=nblk*nthr;

      {
	cl_ulong xalc, xmem;

	clGetDeviceInfo(devID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &xalc, NULL);
	clGetDeviceInfo(devID, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &xmem, NULL);

	while(npho>0){
	  pmax=ntot*npho;
	  pmxo=pmax/OVER;
	  pmax=pmxo*OVER;
	  d.hnum=pmax/hquo;

	  unsigned long aux, mmax=0, mtot=sizeof(datz)+sizeof(dats)+d.gsize*sizeof(DOM);
	  aux=d.hnum*sizeof(hit); if(aux>mmax) mmax=aux; mtot+=aux;
	  aux=pmxo*sizeof(photon); if(aux>mmax) mmax=aux; mtot+=aux;
	  aux=pmax*sizeof(pbuf); if(aux>mmax) mmax=aux; mtot+=aux;

	  if(mmax>xalc || mtot>xmem) npho/=2; else break;
	}
      }

      {
	unsigned int size=d.rsize;
	if(size<ntot) cerr<<"Error: not enough multipliers: only have "<<size<<" (need "<<ntot<<")!"<<endl;
	else d.rsize=ntot;
      }

      unsigned long tot=0, cnt=0;
      cl_int err;

      {
	unsigned long size=sizeof(datz); tot+=size;
	ez = clCreateBuffer(ctx, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, size, &z, &err); checkError(err);
      }

      {
	unsigned long size=d.hnum*sizeof(hit); tot+=size;
	eh = clCreateBuffer(ctx, CL_MEM_READ_WRITE, size, NULL, &err); checkError(err);
      }

      {
	unsigned long size=sizeof(photon);
	size*=pmxo; tot+=size;
	ep = clCreateBuffer(ctx, CL_MEM_READ_WRITE, size, NULL, &err); checkError(err);
      }

      {
	unsigned long size=pmax*sizeof(pbuf); tot+=size;
	bf = clCreateBuffer(ctx, CL_MEM_READ_WRITE, size, NULL, &err); checkError(err);
      }

      {
	unsigned long size=d.gsize*sizeof(DOM); cnt+=size;
	eo = clCreateBuffer(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, size, q.oms, &err); checkError(err);
      }

      {
	unsigned long size=sizeof(dats); tot+=size;
	ed = clCreateBuffer(ctx, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, size, &d, &err); checkError(err);
      }

      cerr<<"Total GPU memory usage: "<<tot<<"  const: "<<cnt<<"  (npho="<<npho<<")"<<endl;

      {
	int n=1;
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &ed));
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &ez));
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &eh));
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &ep));
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &bf));
	checkError(clSetKernelArg(clkernel, n++, sizeof(cl_mem), &eo));
      }
    }

    void fin(){
      fflush(stdout);

      checkError(clReleaseMemObject(ez));
      checkError(clReleaseMemObject(eh));
      checkError(clReleaseMemObject(ep));
      checkError(clReleaseMemObject(bf));
      checkError(clReleaseMemObject(eo));
      checkError(clReleaseMemObject(ed));
    }

    void set(){
      // if(xgpu) device=device;
    }

    void my_wait(){
      if(true){
	cl_event e;
	checkError(clEnqueueMarkerWithWaitList(cq, 0, NULL, &e));

	while(true){
	  cl_int v;
	  checkError(clGetEventInfo(e, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &v, NULL));
	  if(v==CL_COMPLETE) break;
	  usleep(1000);
	}
      }
      else checkError(clFinish(cq));
    }

    void kernel_i(){
      checkError(clEnqueueReadBuffer(cq, ed, CL_FALSE, 0, 2*sizeof(int), &d, 0, NULL, NULL));
      my_wait();

      cl_ulong t1, t2;
      checkError(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t1, NULL));
      checkError(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &t2, NULL));
      deviceTime+=(long long)(t2-t1)/1.e6;

      if(d.ab>0) cerr<<"Error: TOT was a nan or an inf "<<d.ab<<" times!"<<endl;

      if(d.hidx>d.hnum){ cerr<<"Error: data buffer overflow occurred: "<<d.hidx<<">"<<d.hnum<<"!"<<endl; d.hidx=d.hnum; }

      if(d.hidx>0){
	unsigned int size=d.hidx*sizeof(hit);
	checkError(clEnqueueReadBuffer(cq, eh, CL_FALSE, 0, size, &q.hits[xppc::d.hidx], 0, NULL, NULL));
	xppc::d.hidx+=d.hidx;
      }
    }

    void kernel_c(){
      if(old>0) my_wait();
      if(num>0){
	unsigned int size=pk*sizeof(photon);
	checkError(clEnqueueWriteBuffer(cq, ep, CL_FALSE, 0, size, q.pz, 0, NULL, NULL));
      }
    }

    void kernel_f(){
      my_wait();
      if(num>0){
	unsigned int zero=0;
	checkError(clSetKernelArg(clkernel, 0, sizeof(unsigned int), &zero));
	checkError(clEnqueueTask(cq, clkernel, 0, NULL, NULL));

	checkError(clSetKernelArg(clkernel, 0, sizeof(unsigned int), &num));
	checkError(clEnqueueNDRangeKernel(cq, clkernel, 1, NULL, &ntot, &nthr, 0, NULL, &event));
	checkError(clFlush(cq));
      }
    }

    void stop(){
      fprintf(stderr, "Device time: %2.1f [ms]\n", deviceTime);
      if(clkernel) checkError(clReleaseKernel(clkernel));
      if(program) checkError(clReleaseProgram(program));
      if(cq) checkError(clReleaseCommandQueue(cq));
      if(ctx) checkError(clReleaseContext(ctx));
    }

    void ini_f(unsigned int & c, unsigned int t, unsigned int div = 1){
      unsigned int aux=pmax/div;
      d.gini=c, d.gspc=aux, d.gtot=t/div; d.gdiv=div; c+=aux;

      checkError(clEnqueueWriteBuffer(cq, ed, CL_TRUE, 0, 8*sizeof(int), &d, 0, NULL, NULL));

      cerr<<" "<<d.gspc;
    }
  };

  int gcd(int a, int b){
    if(a==0) return b;
    return gcd(b%a, a);
  }

  vector<gpu> gpus;

  void ini(){
    setq();

    d.hnum=0; d.gnum=gpus.size();
    pmax=0, pmxo=0, pn=0, pk=0;

    unsigned int div=1;
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++){
      i->set();
      i->ini(); if(xgpu) sv++;
      d.hnum+=i->d.hnum;
      pmax+=i->pmax;
      if(pmxo==0 || pmxo>i->pmxo) pmxo=i->pmxo;

      if(i->device==0) div=i->pmax;
      else div=gcd(i->pmax, div);
    }

    unsigned int gsum=0; cerr<<"Relative GPU loadings:";
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->ini_f(gsum, pmax, div);
    d.gini=0; d.gspc=gsum; d.gtot=gsum; d.gdiv=div; cerr<<endl;

    {
      q.hits = new hit[d.hnum];
    }

    {
      q.pz = new photon[pmxo];
    }
  }

  void fin(){
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->fin();
    delete q.pz;
    delete [] q.hits;
  }

  void listDevices(){
    cl_device_type dtype=0;
    char * only;

    only=getenv("OCPU");
    if(only!=NULL) if(*only==0 || atoi(only)>0){
	dtype|=CL_DEVICE_TYPE_CPU;
	cerr<<"Running only on CPUs!"<<endl;
    }
    only=getenv("OGPU");
    if(only!=NULL) if(*only==0 || atoi(only)>0){
	dtype|=CL_DEVICE_TYPE_GPU;
	cerr<<"Running only on GPUs!"<<endl;
    }
    only=getenv("OACC");
    if(only!=NULL) if(*only==0 || atoi(only)>0){
	dtype|=CL_DEVICE_TYPE_ACCELERATOR;
	cerr<<"Running only on ACCs!"<<endl;
    }

    if(dtype==0) dtype=CL_DEVICE_TYPE_ALL;

    cl_uint num;
    checkError(clGetPlatformIDs(0, NULL, &num));
    fprintf(stderr, "Found %d platforms:\n", num);

    cl_platform_id ids[num];
    checkError(clGetPlatformIDs(num, ids, NULL));

    for(unsigned int i=0, n=0; i<num; ++i){
      size_t siz;
      cl_platform_id & platform=ids[i];
      fprintf(stderr, "platform %d: ", i);

      {
	checkError(clGetPlatformInfo(platform, CL_PLATFORM_NAME, 0, NULL, &siz));
	char buf[siz+1]; buf[siz]='\0';
	checkError(clGetPlatformInfo(platform, CL_PLATFORM_NAME, siz, &buf, NULL));
	fprintf(stderr, "%s ", buf);
      }

      {
	checkError(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 0, NULL, &siz));
	char buf[siz+1]; buf[siz]='\0';
	checkError(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, siz, &buf, NULL));
	fprintf(stderr, "%s ", buf);
      }

      cl_uint num;
      if(clGetDeviceIDs(platform, dtype, 0, NULL, &num)!=CL_SUCCESS) num=0;
      fprintf(stderr, "device count: %d\n", num);

      if(num>0){
	cl_device_id devices[num];
	checkError(clGetDeviceIDs(platform, dtype, num, devices, NULL));

	for(unsigned int j=0; j<num; j++, n++){
	  cl_device_id & di = devices[j];
	  fprintf(stderr, "  device %d:", j);

	  {
	    cl_uint info;
	    clGetDeviceInfo(di, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &info, NULL);
	    fprintf(stderr, " cu=%d", info);
	  }

	  {
	    size_t info;
	    clGetDeviceInfo(di, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &info, NULL);
	    fprintf(stderr, " gr=%lu", info);
	  }

	  {
	    cl_uint info;
	    clGetDeviceInfo(di, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &info, NULL);
	    fprintf(stderr, " %d MHz", info);
	  }

	  {
	    cl_ulong info;
	    clGetDeviceInfo(di, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &info, NULL);
	    fprintf(stderr, " %llu bytes", (unsigned long long) info);
	  }

	  {
	    cl_ulong info;
	    clGetDeviceInfo(di, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &info, NULL);
	    fprintf(stderr, " cm=%llu", (unsigned long long) info);
	  }

	  {
	    cl_ulong info;
	    clGetDeviceInfo(di, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &info, NULL);
	    fprintf(stderr, " lm=%llu", (unsigned long long) info);
	  }

	  {
	    cl_bool info;
	    clGetDeviceInfo(di, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(cl_bool), &info, NULL);
	    fprintf(stderr, " ecc=%d", info);
	  }

	  all.push_back(make_pair(platform, di));
	  fprintf(stderr, "\n");
	}
      }
    }
    fprintf(stderr, "\n");
  }

  static unsigned int old=0;

  void print();

  void kernel(unsigned int num){
    if(old>0){
      d.hidx=0;
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_i();
      cerr<<"photons: "<<old<<"  hits: "<<d.hidx<<endl;
    }

    {
      unsigned int res=num/d.gspc;
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->num=res*i->d.gspc; res=num-res*d.gspc;
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++){
	unsigned int del=min(res, i->d.gspc);
	i->num+=del; res-=del;
      }
    }

    {
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_c();
    }

    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_f();

    if(old>0) print();

    old=num;
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->old=i->num;
  }

  void start(){
  }

  void stop(){
    fprintf(stderr, "\n");
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->stop();
  }

  void choose(int device){
    listDevices();
    int deviceCount=all.size();

    if(device<0){
      if(deviceCount==0){ cerr<<"Could not find compatible devices!"<<endl; exit(2); }
      for(int device=0; device<deviceCount; ++device){
	gpus.push_back(gpu(device));
	if(gpus.back().npho<=0) gpus.pop_back();
      }
    }
    else{
      if(device>=deviceCount){ cerr<<"Device #"<<device<<" is not available!"<<endl; exit(2); }
      sv+=device;
      gpus.push_back(gpu(device));
      if(gpus.back().npho<=0) gpus.pop_back();
    }
    if(gpus.size()<=0){
      cerr<<"No active GPU(s) selected!"<<endl;
      exit(2);
    }
    xgpu=gpus.size()>1;
  }

#include "f2k.cxx"
}

#ifndef XLIB
using namespace xppc;

float zshift(cl_float4 r){
  float edh=d.dh;
  return zshift(&d, r, 0, &edh, &z);
}

int main(int arg_c, char *arg_a[]){
  cerr<<"sizes: DOM="<<sizeof(DOM)<<" hit="<<sizeof(hit)<<" pbuf="<<sizeof(pbuf)<<" photon="<<sizeof(photon)<<endl;
  start();
  if(arg_c<=1){
    listDevices();
    fprintf(stderr, "Use: %s [device] (f2k muons)\n"
	    "     %s [str] [om] [num] [device] (flasher)\n", arg_a[0], arg_a[0]);
  }
  else if(0==strcmp(arg_a[1], "-")){
    initialize();
    ices & w = z.w[WNUM/2];
    cerr<<"For wavelength="<<q.wvs[w.wvl]<<" [nm]  np="<<(1/w.coschr)<<"  cm="<<1/w.ocm<<" [m/ns]"<<endl;
    cl_float4 r;
    r.w=0;
    if(arg_c==4){
      r.x=atof(arg_a[2]);
      r.y=atof(arg_a[3]);
    }
    else r.x=0, r.y=0;
    for(int i=0; i<d.size; i++){
      float z=d.hmin+d.dh*i;
      r.z=z; for(int j=0; j<10; j++) r.z=z+zshift(r); z=r.z;
      cout<<z<<" "<<w.z[i].abs<<" "<<w.z[i].sca*(1-d.g)<<" "<<d.az[i].ra*d.sum<<endl;
    }
  }
  else if(0==strcmp(arg_a[1], "_")){
    initialize();
    cl_float4 r;
    r.w=0;
    for(r.x=-750.f; r.x<751.f; r.x+=3.f) for(r.y=-750.f; r.y<751.f; r.y+=3.f) for(float z=-750.f; z<751.f; z+=6.f){
	  r.z=z; for(int j=0; j<10; j++) r.z=z+zshift(r);
	  cout<<z<<" "<<r.x<<" "<<r.y<<" "<<(r.z-z)<<endl;
	}
  }
  else if(arg_c<=2){
    int device=0;
    if(arg_c>1) device=atoi(arg_a[1]);
    initialize();
    choose(device);
    fprintf(stderr, "Processing f2k muons from stdin on device %d\n", device);
    f2k();
  }
  else{
    int str=0, dom=0, device=0, itr=0;
    unsigned long long num=1000000ULL;

    if(arg_c>1) str=atoi(arg_a[1]);
    if(arg_c>2) dom=atoi(arg_a[2]);
    if(arg_c>3){
      num=(unsigned long long) atof(arg_a[3]);
      char * sub = strchr(arg_a[3], '*');
      if(sub!=NULL) itr=(int) atof(++sub);
    }
    if(arg_c>4) device=atoi(arg_a[4]);
    initialize();
    choose(device);
    fprintf(stderr, "Running flasher simulation on device %d\n", device);
    flasher(str, dom, num, itr);
  }

  stop();
}
#endif
