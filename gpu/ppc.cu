#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sys/time.h>

#ifndef __CUDACC__
#define XCPU
#endif

#ifdef XCPU
#include <cmath>
#include <cstring>
#endif

using namespace std;

namespace xppc{
#include "ini.cxx"
#include "pro.cu"

  void initialize(float enh = 1.f){ m.set(); q.eff*=enh; }

  unsigned int pmax, pmxo, pn, pk, hquo;

  void setq(){
    char * HQUO=getenv("HQUO");
    hquo=HQUO==NULL?1:atoi(HQUO);
    cerr<<"HQUO(photons/max number of hits)="<<hquo<<endl;
  }

#ifdef XCPU
  dats *e;  // pointer to a copy of "d" on device
  int nblk, nthr, ntot;

  void ini(){
    setq();
    rs_ini();
    pn=0; pk=0;

    ntot=nblk*nthr;
    pmax=ntot*NPHO;
    pmxo=pmax/OVER;
    pmax=pmxo*OVER;
    d.hnum=pmax/hquo;

    d.gdev=0; d.gnum=1;
    d.gini=0; d.gspc=pmax; d.gtot=pmax; d.gdiv=1;

    {
      d.hits = q.hits = new hit[d.hnum];
      d.pz = q.pz = new photon[pmxo];
      d.bf = new pbuf[pmax];
    }

    {
      d.z=&z; d.oms=q.oms; e=&d;
    }

    {
      unsigned int size=d.rsize, need=seed+1;
      if(size<need) cerr<<"Error: not enough multipliers: asked for "<<seed<<"-th out of "<<size<<"!"<<endl;
    }
  }

  void fin(){
    delete d.pz;
    delete d.hits;
    delete d.bf;
  }
#else
  bool xgpu=false;

  void checkError(cudaError result){
    if(result!=cudaSuccess){
      cerr<<"CUDA Error: "<<cudaGetErrorString(result)<<endl;
      exit(2);
    }
  }

  struct gpu{
    dats d;
    dats *e;  // pointer to a copy of "d" on device

    int device;
    int nblk, nthr, ntot;
    unsigned int npho, pmax, pmxo;

    float dt, deviceTime, threadMin, threadMax;
    cudaDeviceProp prop;
    cudaStream_t stream;
    cudaEvent_t evt1, evt2;

    unsigned int old, num;

    gpu(int device) : deviceTime(0), threadMin(0), threadMax(0), old(0), npho(NPHO){
      this->device=device;

      {
	ostringstream o; o<<"NPHO_"<<device;
	char * nph=getenv(o.str().c_str());
	if(nph==NULL) nph=getenv("NPHO");
	if(nph!=NULL) if(*nph!=0){
	    npho=atoi(nph);
	    cerr<<"Setting NPHO="<<npho<<endl;
	    if(npho<=0){
	      cerr<<"Not using device # "<<device<<"!"<<endl;
	      return;
	    }
	  }
      }

      checkError(cudaSetDevice(device));
      checkError(cudaGetDeviceProperties(&prop, device));

#if CUDART_VERSION >= 3000
      checkError(cudaFuncSetCacheConfig(propagate, cudaFuncCachePreferL1));
#endif

      cudaFuncAttributes attr;
      checkError(cudaFuncGetAttributes (&attr, propagate));

      nblk=prop.multiProcessorCount, nthr=attr.maxThreadsPerBlock;
      cerr<<"Running on "<<nblk<<" MPs x "<<nthr<<" threads"<<endl;

      fprintf(stderr, "Kernel uses: l=%lu r=%d s=%lu c=%lu\n", (unsigned long)attr.localSizeBytes,
	      attr.numRegs, (unsigned long)attr.sharedSizeBytes, (unsigned long)attr.constSizeBytes);
    }

    void ini(){
      rs_ini();
      d=xppc::d;

      d.gdev=device;

      {
	d.blockIdx=-1, d.gridDim=nblk;
	ostringstream o; o<<"BADMP_"<<device;
	char * bmp=getenv(o.str().c_str());
	if(bmp==NULL) bmp=getenv("BADMP");
	if(bmp!=NULL) if(*bmp!=0){
	    d.blockIdx=atoi(bmp), d.gridDim--;
	    cerr<<"Not using MP #"<<d.blockIdx<<endl;
	  }
      }

      ntot=nblk*nthr;

      {
	unsigned long xmem=prop.totalGlobalMem;

	while(npho>0){
	  pmax=ntot*npho;
	  pmxo=pmax/OVER;
	  pmax=pmxo*OVER;
	  d.hnum=pmax/hquo;

	  unsigned long mtot=sizeof(datz)+sizeof(dats)+d.gsize*sizeof(DOM);
	  mtot+=d.hnum*sizeof(hit);
	  mtot+=pmxo*sizeof(photon);
	  mtot+=pmax*sizeof(pbuf);

	  if(mtot>xmem) npho/=2; else break;
	}
      }

      {
	checkError(cudaStreamCreate(&stream));
	checkError(cudaEventCreateWithFlags(&evt1, cudaEventBlockingSync));
	checkError(cudaEventCreateWithFlags(&evt2, cudaEventBlockingSync));
      }

      {
	unsigned int size=d.rsize;
	if(size<ntot) cerr<<"Error: not enough multipliers: only have "<<size<<" (need "<<ntot<<")!"<<endl;
	else d.rsize=ntot;
      }

      unsigned long tot=0, cnt=0;

      {
	unsigned long size=sizeof(datz); tot+=size;
	checkError(cudaMalloc((void**) &d.z, size));
	checkError(cudaMemcpy(d.z, &z, size, cudaMemcpyHostToDevice));
      }

      {
	unsigned long size=d.hnum*sizeof(hit); tot+=size;
	checkError(cudaMalloc((void**) &d.hits, size));
      }

      {
	unsigned long size=pmxo*sizeof(photon); tot+=size;
	checkError(cudaMalloc((void**) &d.pz, size));
      }

      {
	unsigned long size=pmax*sizeof(pbuf); tot+=size;
	checkError(cudaMalloc((void**) &d.bf, size));
      }

      {
	unsigned long size=d.gsize*sizeof(DOM); cnt+=size;
	checkError(cudaMalloc((void**) &d.oms, size));
	checkError(cudaMemcpy(d.oms, &q.oms, size, cudaMemcpyHostToDevice));
	// checkError(cudaMemcpyToSymbol(oms, q.oms, size));
      }

      {
	unsigned long size=sizeof(dats); tot+=size;
	checkError(cudaMalloc((void**) &e, size));
	checkError(cudaMemcpy(e, &d, size, cudaMemcpyHostToDevice));
      }

      cerr<<"Total GPU memory usage: "<<tot<<"  const: "<<cnt<<"  (npho="<<npho<<")"<<endl;
    }

    void fin(){
      checkError(cudaFree(d.z));
      checkError(cudaFree(d.hits));
      checkError(cudaFree(d.pz));
      checkError(cudaFree(d.bf));
      checkError(cudaFree(d.oms));
      checkError(cudaFree(e));
      checkError(cudaEventDestroy(evt1));
      checkError(cudaEventDestroy(evt2));
      checkError(cudaStreamDestroy(stream));
    }

    void set(){
      if(xgpu) checkError(cudaSetDevice(device));
    }

    void kernel_i(){
      {
	checkError(cudaStreamSynchronize(stream));
	checkError(cudaMemcpy(&d, e, 7*sizeof(int), cudaMemcpyDeviceToHost));

	checkError(cudaEventElapsedTime(&dt, evt1, evt2)); deviceTime+=dt;

	if(d.ab>0){
	  cerr<<"Error: TOT was a nan or an inf "<<d.ab<<" times! Bad GPU "<<device<<" MP";
	  for(int i=0; i<min(d.ab, 4); i++) cerr<<" #"<<d.bmp[i]; cerr<<endl;
	}
	if(d.mp!=d.gridDim){ cerr<<"Error: did not encounter MP #"<<d.blockIdx<<endl; exit(4); }
	if(threadMax!=-1){
	  if((unsigned long long)(dt*prop.clockRate)<0x100000000ULL){
	    threadMin+=d.tn/(float)prop.clockRate;
	    threadMax+=d.tx/(float)prop.clockRate;
	  }
	  else threadMin=-1, threadMax=-1;
	}

	if(d.hidx>d.hnum){ cerr<<"Error: data buffer overflow occurred: "<<d.hidx<<">"<<d.hnum<<"!"<<endl; d.hidx=d.hnum; }
      }

      {
	unsigned int size=d.hidx*sizeof(hit);
	checkError(cudaMemcpyAsync(&q.hits[xppc::d.hidx], d.hits, size, cudaMemcpyDeviceToHost, stream));
	xppc::d.hidx+=d.hidx;
      }
    }

    void kernel_c(){
      if(old>0) checkError(cudaStreamSynchronize(stream));
      if(num>0){
	unsigned int size=pk*sizeof(photon);
	checkError(cudaMemcpyAsync(d.pz, q.pz, size, cudaMemcpyHostToDevice, stream));
      }
    }

    void kernel_f(){
      checkError(cudaStreamSynchronize(stream));
      if(num>0){
	checkError(cudaEventRecord(evt1, stream));

	propagate<<< 1, 1, 0, stream >>>(e, 0);
	checkError(cudaGetLastError());

	propagate<<< nblk, nthr, 0, stream >>>(e, num);
	checkError(cudaGetLastError());

	checkError(cudaEventRecord(evt2, stream));
      }
    }

    void stop(){
      fprintf(stderr, "Device time: %2.1f (in-kernel: %2.1f...%2.1f) [ms]\n", deviceTime, threadMin, threadMax);
      checkError(cudaDeviceReset());
    }

    void ini_f(unsigned int & c, unsigned int t, unsigned int div = 1){
      unsigned int aux=pmax/div;
      d.gini=c, d.gspc=aux, d.gtot=t/div; d.gdiv=div; c+=aux;

      checkError(cudaMemcpy(e, &d, 14*sizeof(int), cudaMemcpyHostToDevice));

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
    pmax=0, pmxo=0, pn=0; pk=0;

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
      unsigned long size=d.hnum*sizeof(hit);
      checkError(cudaMallocHost((void**) &q.hits, size));
    }

    {
      unsigned long size=pmxo*sizeof(photon);
      checkError(cudaMallocHost((void**) &q.pz, size));
    }
  }

  void fin(){
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->fin();
    checkError(cudaFreeHost(q.hits));
    checkError(cudaFreeHost(q.pz));
  }

  void listDevices(){
    int deviceCount, driver, runtime;
    cudaGetDeviceCount(&deviceCount);
    cudaDriverGetVersion(&driver);
    cudaRuntimeGetVersion(&runtime);
    fprintf(stderr, "Found %d devices, driver %d, runtime %d\n", deviceCount, driver, runtime);
    for(int device=0; device<deviceCount; ++device){
      cudaDeviceProp prop; cudaGetDeviceProperties(&prop, device);
      fprintf(stderr, "%d(%d.%d): %s %g GHz G(%lu) S(%lu) C(%lu) R(%d) W(%d)\n"
	      "\tl%d o%d c%d h%d i%d m%d a%lu M(%lu) T(%d: %d,%d,%d) G(%d,%d,%d)\n",
	      device, prop.major, prop.minor, prop.name, prop.clockRate/1.e6,
	      (unsigned long)prop.totalGlobalMem, (unsigned long)prop.sharedMemPerBlock,
	      (unsigned long)prop.totalConstMem, prop.regsPerBlock, prop.warpSize,
	      prop.kernelExecTimeoutEnabled, prop.deviceOverlap, prop.computeMode,
	      prop.canMapHostMemory, prop.integrated, prop.multiProcessorCount,
	      (unsigned long)prop.textureAlignment,
	      (unsigned long)prop.memPitch, prop.maxThreadsPerBlock,
	      prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2],
	      prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    }
    fprintf(stderr, "\n");
  }

  static unsigned int old=0;
#endif

  void print();

  void kernel(unsigned int num){
#ifdef XCPU
    unsigned int & old = num;
#endif
    if(old>0){
      d.hidx=0;
#ifdef XCPU
      for(d.blockIdx=0, d.gridDim=nblk, blockDim.x=nthr; d.blockIdx<d.gridDim; d.blockIdx++)
	for(threadIdx.x=0; threadIdx.x<blockDim.x; threadIdx.x++) propagate(e, num);

      if(d.hidx>d.hnum){ cerr<<"Error: data buffer overflow occurred: "<<d.hidx<<">"<<d.hnum<<"!"<<endl; d.hidx=d.hnum; }
#else
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_i();
#endif
      cerr<<"photons: "<<old<<"  hits: "<<d.hidx<<endl;
    }

#ifndef XCPU
    {
      unsigned int res=num/d.gspc;
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->num=res*i->d.gspc; res=num-res*d.gspc;
      for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++){
	unsigned int del=min(res, i->d.gspc);
	i->num+=del; res-=del;
      }
    }

    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_c();

    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->kernel_f();
#endif

    if(old>0) print();
#ifndef XCPU
    old=num;
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->old=i->num;
#endif
  }

#ifdef XCPU
  void start(){}
  void stop(){}
  void choose(int device){
    sv+=device;
    seed=device;
    nblk=NBLK, nthr=NTHR;
  }
  void listDevices(){}
#else
  void start(){
    cudaSetDeviceFlags(cudaDeviceBlockingSync);
  }

  void stop(){
    fprintf(stderr, "\n");
    for(vector<gpu>::iterator i=gpus.begin(); i!=gpus.end(); i++) i->set(), i->stop();
  }

  void choose(int device){
    if(device<0){
      int deviceCount; cudaGetDeviceCount(&deviceCount);
      for(int device=0; device<deviceCount; ++device){
	gpus.push_back(gpu(device));
	if(gpus.back().npho<=0) gpus.pop_back();
      }
    }
    else{
      sv+=device;
      gpus.push_back(gpu(device));
      if(gpus.back().npho<=0) gpus.pop_back();
    }
    if(gpus.size()<=0){
      cerr<<"No active GPU(s) selected!"<<endl;
      exit(5);
    }
    xgpu=gpus.size()>1;
  }
#endif

#include "f2k.cxx"
}

#ifndef XLIB
using namespace xppc;

float zshift(float4 r){
  float edh=d.dh;
  return zshift(d, r, 0, edh);
}

int main(int arg_c, char *arg_a[]){
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
    float4 r;
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
    float4 r;
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
