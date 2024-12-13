#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sys/time.h>

#include <cmath>
#include <cstring>

#include <atomic>
#include <iterator>
#include <execution>
#include <thrust/iterator/counting_iterator.h>

#ifdef USE_I3_LOGGING
#include "icetray/I3Logging.h"
#else
#define log_info_stream(msg) \
  do { std::cerr << msg << std::endl; } while (0)
#endif

using namespace std;

namespace xppc{
  void sincosf(float x, float * s, float * c){ *s=sinf(x), *c=cosf(x); }

#include "ini.cxx"
#include "pro.cxx"

  void initialize(float enh = 1.f){ m.set(); q.eff*=enh; }

  unsigned int pmax, pmxo, pn, pk, hquo;

  void setq(){
    char * HQUO=getenv("HQUO");
    hquo=HQUO==NULL?1:atoi(HQUO);
    cerr<<"HQUO(photons/max number of hits)="<<hquo<<endl;
  }

  dats *e;  // pointer to a copy of "d" on device

  void ini(){
    setq();
    rs_ini();
    pn=0; pk=0;

    dats & d = * ed;

    int ntot;
    char * NTHR=getenv("NTHR");
    ntot=NTHR==NULL?NTOT:atoi(NTHR);
    cerr<<"NTHR(threads)="<<ntot<<endl;
    d.ntot=ntot;

    pmax=ntot*NPHO;
    pmxo=pmax/OVER;
    pmax=pmxo*OVER;
    d.hnum=pmax/hquo;

    {
      q.hits = new hit[d.hnum];
      q.pz = new photon[pmxo];
      q.bf = new pbuf[pmax];
    }

    {
      oms = new DOM[MAXGEO];
      memcpy(oms, q.oms, MAXGEO*sizeof(DOM));
    }

    {
      unsigned int size=d.rsize;
      if(size<ntot) cerr<<"Error: asked for more threads ("<<ntot<<") than multipliers ("<<size<<") !"<<endl;
    }
  }

  void fin(){
    delete q.pz;
    delete q.hits;
    delete q.bf;
    delete oms;
  }

#ifdef XLIB
  size_t getMaxBunchSize(){ return pmxo; }
  size_t getWorkgroupSize(){ return ed->ntot; }
  float getTotalDeviceTime(){ return 0.f; }
#endif

  void print();

  void kernel(unsigned int num){
    if(num>0){
      dats & d = * ed;
      d.hidx=0;

      propagate(num, ed, ez, q.hits, q.pz, q.bf, oms);

      if(d.hidx>d.hnum){ cerr<<"Error: data buffer overflow occurred: "<<d.hidx<<">"<<d.hnum<<"!"<<endl; d.hidx=d.hnum; }

      log_info_stream("photons: "<<num<<"  hits: "<<d.hidx);
      print();
    }
  }

  void start(){}
  void stop(){}
  void choose(int device){
    sv+=device;
  }
  void listDevices(){}

#include "f2k.cxx"
}

#ifndef XLIB
using namespace xppc;

float zdh;

float zshift(float4 r){
  dats & d = * ed;
  zdh=d.dh;
  return zshift(d, r, zdh, *ez);
}

int main(int arg_c, char *arg_a[]){
  dats & d = * ed;
  datz & z = * ez;

  start();
  if(arg_c<=1){
    listDevices();
    fprintf(stderr, "Use: %s [device] (f2k muons)\n"
	    "     %s [str] [om] [num] [device] (flasher)\n", arg_a[0], arg_a[0]);
  }
  else if(0==strcmp(arg_a[1], "-")){
    initialize();
    ices & w = z.w[WNUM/2];
    cerr<<"For wavelength="<<q.wvs[w.wvl].w<<" [nm]  np="<<(1/w.coschr)<<"  cm="<<1/w.ocm<<" [m/ns]"<<endl;
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
  else if(0==strcmp(arg_a[1], "=")){
    initialize();
    ices & w = z.w[WNUM/2];
    cerr<<"For wavelength="<<q.wvs[w.wvl].w<<" [nm]  np="<<(1/w.coschr)<<"  cm="<<1/w.ocm<<" [m/ns]"<<endl;
    float4 r;
    r.w=0;
    string in;
    while(getline(cin, in)){
      if(3==sscanf(in.c_str(), "%f %f %f", &r.x, &r.y, &r.z)){
	float dz=zshift(r);
	cout<<in<<" "<<dz<<" "<<zdh<<endl;
      }
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
