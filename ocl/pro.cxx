#ifdef SHRT
#define STRINGIFY(A) #A
#define XTRINGIFY(A) STRINGIFY(A)
string kernel_source = XTRINGIFY((
#elif __APPLE_CC__
typedef union{
  cl_float  CL_ALIGNED(16) s[4];
  __extension__ struct{ cl_float x, y, z, w; };
} cl_float4;
#endif

typedef struct{
  cl_float R, F;
  cl_float r[3];
} DOM;

typedef struct{
  cl_uint i, n, z;
  cl_float t;
  cl_float pth, pph, dth, dph;
} hit;

typedef struct{
  cl_float4 r;    // location, time
  cl_float4 n;    // direction
  cl_uint q;      // wavelength bin
  cl_uint i;      // segment number
  cl_int fla, ofla;
} pbuf;

typedef struct{
  cl_float4 r;    // location, time
  cl_float4 n;    // direction, track length
  cl_uint q;      // track segment
  cl_uint num;    // number of photons in this bunch
  cl_int type;    // source type
  cl_float f;     // fraction of light from muon alone (without cascades)
  union{
    struct{
      cl_float a, b;  // longitudinal development parametrization coefficients
      cl_float beta;  // velocity of incident particle
      cl_float tau;   // luminescence decay time
    };
    struct{
      cl_float ka, up; // 2d-gaussian rms and zenith of cone
      cl_float fldr;   // horizontal direction of the flasher led #1
      cl_short fla, ofla;
    };
    cl_int4 c;  // for quick copy
  };
} photon;

typedef struct{
  cl_int wvl;               // wavelength number of this block
  cl_float ocm;             // 1 / speed of light in medium
  cl_float coschr, sinchr;  // cos and sin of the cherenkov angle
  struct{
    cl_float abs;           // absorption
    cl_float sca;           // scattering
  } z [MAXLYS];
} ices;

typedef struct{
  cl_float k1;
  cl_float k2;
  cl_float ra;
  cl_float rb;
} aniz;

typedef struct{
  cl_short n, max;
  cl_float x, y, r;
  cl_float h, d;
  cl_float dl, dh;
} line;

typedef struct{
  ices w[WNUM];
  cl_float2 lp[LMAX][LYRS];
  cl_uint rm[MAXRND];
  cl_ulong rs[MAXRND];
} datz;

typedef struct{
  cl_uint hidx;
  cl_uint ab;    // if TOT was abnormal
  cl_uint gdev;  // number of this GPU
  cl_uint gnum;  // number of all GPUs

  cl_uint gini, gspc, gtot, gdiv;

  cl_float rx;
  cl_float hifl;

  cl_uint hnum;  // size of hits buffer
  cl_int size;   // size of kurt table
  cl_int rsize;  // count of multipliers
  cl_int gsize;  // count of initialized OMs

  cl_short tmod; // tilt model: 1: 1d, 2: 2d
  cl_short vthk; // 0: uniform layers, 1: variable thickness
  cl_float dh, rdh, hmin; // step, 1/step, and min depth

  cl_float ocv;  // 1 / speed of light in vacuum
  cl_float sf;   // scattering function: 0=HG; 1=SAM
  cl_float g, gr; // g=<cos(scattering angle)> and gr=(1-g)/(1+g)

  cl_float xR;   // DOM oversize scaling factor
  cl_float SF, G, GR; // hole ice sf, g, gr
  cl_float hr, hr2, hs, ha; // hole ice radius, radius^2, effective scattering and absorption coefficients

  cl_float azx, azy;  // ice anisotropy direction

  cl_int cn[2];
  cl_float cl[2], crst[2];

  cl_float cb[2][2];

  cl_int lnum, lpts, l0;
  cl_float lmin, lrdz, r0;
  cl_float lnx, lny;
  cl_float lr[LMAX];

  cl_float mmin[2], mstp[2];
  cl_int mnum[2], mcut[2];

  cl_float k1, k2, kz, fr; // ice absorption anisotropy parameters
  aniz az[MAXLYS];

  cl_ushort ls[NSTR];
  cl_ushort is[CX][CY];
  cl_char mcol[CTX][CTY];

  cl_float sum, bfr[12];

  line sc[NSTR];
} dats;

#ifdef SHRT
#define sin native_sin
#define cos native_cos
#define pow native_powr
#define exp native_exp
#define log native_log
#define exp2 native_exp2
#define sqrt native_sqrt
#define rsqrt native_rsqrt

#ifndef CL_VERSION_1_2
#define clamp(x,y,z) min(max(x,y),z)
#endif

float xrnd(uint4 * s){
  ulong sda = (*s).z * (ulong) (*s).x + (*s).y;
  (*s).x = sda; (*s).y = sda >> 32;
  uint tmp = (*s).x >> 9;
  return 2.0f-as_float(tmp|0x3f800000);
}

float mrnd(float k, uint4 * s){  // gamma distribution
  float x;
  if(k<1){  // Weibull algorithm
    float c=1/k;
    float d=(1-k)*pow(k, 1/(c-1));
    float z, e;
    do{
      z=-log(xrnd(s));
      e=-log(xrnd(s));
      x=pow(z, c);
    } while(z+e<d+x);
  }
  else{  // Cheng's algorithm
    float b=k-log(4.0f);
    float l=sqrt(2*k-1);
    float c=1+log(4.5f);
    float u, v, y, z, r;
    do{
      u=xrnd(s); v=xrnd(s);
      y=-log(1/v-1)/l;
      x=k*exp(y);
      z=u*v*v;
      r=b+(k+l)*y-x;
    } while(r<4.5f*z-c && r<log(z));
  }
  return x;
}

float grnd(uint4 * s){  // gaussian distribution
  return sqrt(-2*log(xrnd(s)))*sin(2*FPI*xrnd(s));
}

float2 grnd2(uint4 * s){  // gaussian distribution
  float r=sqrt(-2*log(xrnd(s))), xi=2*FPI*xrnd(s);
  float2 p = r*(float2)(cos(xi), sin(xi));
  return p;
}

float3 mx_normalize(float3 n){
  n*=rsqrt(n.x*n.x+n.y*n.y+n.z*n.z);
  return n;
}

float4 my_normalize(float4 n){
  n.xyz*=rsqrt(n.x*n.x+n.y*n.y+n.z*n.z);
  return n;
}

float4 turn(float cs, float si, float4 n, uint4 * s){
  float4 r = (float4)(n.xyz*n.xyz, 0);

  float4 p1 = (float4)(0);
  float4 p2 = (float4)(0);

  if(r.y>r.z){
    if(r.y>r.x){
      p1 = (float4)(n.y, -n.x, 0, 0);
      p2 = (float4)(0, -n.z, n.y, 0);
    }
    else{
      p1 = (float4)(-n.y, n.x, 0, 0);
      p2 = (float4)(-n.z, 0, n.x, 0);
    }
  }
  else{
    if(r.z>r.x){
      p1 = (float4)(n.z, 0, -n.x, 0);
      p2 = (float4)(0, n.z, -n.y, 0);
    }
    else{
      p1 = (float4)(-n.y, n.x, 0, 0);
      p2 = (float4)(-n.z, 0, n.x, 0);
    }
  }

  p1 = my_normalize(p1);
  p2 = my_normalize(p2);

  r = p1-p2; p2 += p1;
  p1 = my_normalize(r);
  p2 = my_normalize(p2);

  float xi = 2*FPI*xrnd(s);
  float2 p = (float2)(cos(xi), sin(xi));

  return my_normalize((float4)(cs*n.xyz+si*(p.x*p1.xyz+p.y*p2.xyz), n.w));
}

#else
#define __local
#define __global
#define convert_int_sat_rtn (int)floor
#define clamp(x,y,z) min(max(x,y),z)
#define int2 cl_int2
#define int3 cl_int3
#define float2 cl_float2
#define float3 cl_float3
#define float4 cl_float4
#endif

float square(float x){
  return x*x;
}

int findcol(__local dats * d, int x, int y){
  return x>=0 && x<d->mnum[0] && y>=0 && y<d->mnum[1] ? d->mcol[x][y] : -1;
}

float midp(float2 q, float z){
  return q.x+q.y*z;
}

float interpolate(__local dats * d, int j, float nr, float z, __global datz * p){
  int i=j-1, k=clamp(convert_int_sat_rtn(z), 0, d->lpts);
  return ( midp(p->lp[j][k], z-k)*(nr-d->lr[i]) +
	   midp(p->lp[i][k], z-k)*(d->lr[j]-nr) )/(d->lr[j]-d->lr[i]);
}

float interpolatE(__local dats * d, int3 c, float3 w, float z, __global datz * p){
  int k=clamp(convert_int_sat_rtn(z), 0, d->lpts);
  return
    ( c.x<0?0: midp(p->lp[c.x][k], z-k)*w.x )+
    ( c.y<0?0: midp(p->lp[c.y][k], z-k)*w.y )+
    ( c.z<0?0: midp(p->lp[c.z][k], z-k)*w.z );
}

float zshift(__local dats * d, float4 r, int sign, float * dz, __global datz * p){
  if(d->tmod==1){
    float nr=d->lnx*r.x+d->lny*r.y-d->r0;
    int j=1; for(; j<LMAX; j++) if(nr<d->lr[j] || j==d->lnum-1) break;

    float z1=(r.z-d->lmin)*d->lrdz;
    float r1=interpolate(d, j, nr, z1, p);
    if(sign!=0){
      float z2=z1+sign**dz*d->lrdz;
      float r2=interpolate(d, j, nr, z2, p);
      *dz /= 1-sign*(r2-r1)/ *dz;
    }
    return r1;
  }
  else if(d->tmod==2){
    float2 q; // ctr(d, r, q);
    q.x=d->cb[0][0]*r.x+d->cb[1][0]*r.y;
    q.y=d->cb[0][1]*r.x+d->cb[1][1]*r.y;

    q.x=(q.x-d->mmin[0])/d->mstp[0];
    q.y=(q.y-d->mmin[1])/d->mstp[1];

    { // extrapolate
      float x=q.x, y=q.y, xy;
      int n0=d->mnum[0]-1, n1=d->mnum[1]-1, c0=d->mcut[0]-1, c1=d->mcut[1]-1;

      if(y<0) if(x-y>0 && x<=c0){ x=(x-y)*c0/(c0-y); y=0; } // 1
      if(x<0) if(y-x>=0 && y<c1){ y=(y-x)*c1/(c1-x); x=0; } // 2

      if(y>n1) if(x>=n1-c1 && y-x>n1-n0){ x=n1-c1+(n0-(n1-c1))*(x-(n1-c1))/(y-(n1-n0)-(n1-c1)); y=n1; } // 5
      if(x>n0) if(y>n0-c0 && y-x<=n1-n0){ y=n0-c0+(n1-(n0-c0))*(y-(n0-c0))/(x-(n0-n1)-(n0-c0)); x=n0; } // 6

      if(y-x>c1) if(y>=c1 && x<n1-c1){ xy=y-x-c1; xy=c1+(x+y-c1+xy)*(n1-c1)/(n1-c1+xy); x=(xy-c1)/2; y=(xy+c1)/2; } // 3
      if(x-y>c0) if(x>c0 && y<=n0-c0){ xy=x-y-c0; xy=c0+(x+y-c0+xy)*(n0-c0)/(n0-c0+xy); x=(xy+c0)/2; y=(xy-c0)/2; } // 4

      q.x=x, q.y=y;
    }

    int2 qn;
    qn.x=clamp(convert_int_sat_rtn(q.x), 0, d->mnum[0]);
    qn.y=clamp(convert_int_sat_rtn(q.y), 0, d->mnum[1]);

    float2 qr;
    qr.x=q.x-qn.x;
    qr.y=q.y-qn.y;

    float qrx, qry;
    int qnx, qny;

    if(qr.x>qr.y){
      qrx=qr.x;
      qry=qr.y;
      qnx=qn.x+1;
      qny=qn.y;
    }
    else{
      qrx=qr.y;
      qry=qr.x;
      qnx=qn.x;
      qny=qn.y+1;
    }

    int3 c;
    c.x=findcol(d, qn.x, qn.y);
    c.y=findcol(d, qnx, qny);
    c.z=findcol(d, qn.x+1, qn.y+1);

    float3 w;
    w.x=1-qrx, w.y=qrx-qry, w.z=qry;

    float z1=(r.z-d->lmin)*d->lrdz;
    float r1=interpolatE(d, c, w, z1, p);
    if(sign!=0){
      float z2=z1+sign**dz*d->lrdz;
      float r2=interpolatE(d, c, w, z2, p);
      *dz /= 1-sign*(r2-r1)/ *dz;
    }

    return r1;
  }
  else return 0;
}

#ifndef SHRT
#undef __local
#undef __global
#undef convert_int_sat_rtn
#undef clamp
#undef int2
#undef int3
#undef float2
#undef float3
#undef float4
#else

float2 ctr(__local dats * d, float2 r){
  return (float2)(d->cb[0][0]*r.x+d->cb[1][0]*r.y, d->cb[0][1]*r.x+d->cb[1][1]*r.y);
}

__kernel void propagate(__private uint num,
			__global dats * ed,
			__global datz * ez,
			__global hit * eh,
			__global photon * ep,
			__global pbuf * bf,
			__global DOM * oms){

  if(num==0){
    if(get_global_size(0)>1) return;

    ed->hidx=0;
    ed->ab=0;

    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    return;
  }

  __local dats e;
  __global ices * w;

  {
    event_t ev=async_work_group_copy((__local char *) &e, (__global char *) ed, sizeof(e), 0);
    wait_group_events(1, &ev);
    if(get_local_id(0)==0) e.hidx=get_global_size(0)+get_group_id(0);
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  const unsigned int idx=get_local_id(0)*get_num_groups(0)+get_group_id(0);

  int niw=0, old=idx%e.rsize;
  uint4 s = (uint4)(ez->rs[old], ez->rs[old] >> 32, ez->rm[old], old);

  float4 r=(float4)(0);
  float4 n=(float4)(0);

  float TOT=0, SCA;

  for(unsigned int i=idx, pj=-1U, pn=0; i<num; i+=get_global_size(0)){
    while(e.gini+i%e.gspc+(i/e.gspc)*e.gtot>=pn) pn+=ep[++pj].num;

    unsigned int j=min(convert_int_sat_rtn(WNUM*xrnd(&s)), WNUM-1);
    w=&ez->w[j];
    n.w=w->ocm;

    photon p=ep[pj];
    r=p.r, n.xyz=p.n.xyz;
    float l=p.n.w; niw=p.q;
    int fla, ofla;

    if(p.type>0){
      fla=p.fla, ofla=p.ofla;

      if(p.type<=4){
	float xi=xrnd(&s);
	if(p.fldr<0) xi*=2*FPI;
	else{
	  int r=convert_int_sat_rtn(p.fldr/360)+1;
	  int s=convert_int_sat_rtn(xi*r);
	  xi=radians(p.fldr+s*360/r);
	}
	n.x=cos(xi), n.y=sin(xi);

	float FLZ=0.07735f, FLR=0.119f-0.008f, FLB=0.03396f;
	if(p.ka>0) r=(float4)(FLR*n.x, FLR*n.y, FLZ, 0);

	float np=cos(p.up); n.z=sin(p.up);
	n.x*=np; n.y*=np;

	if(p.ka>999.f){
	  float pf=xrnd(&s);

	  float   s1=0.0191329f, s2=0.0686944f, C1=0.583583f, C2=0.967762f;
	  switch(p.type){
	  case 1: s1=0.0196242f, s2=0.0750278f, C1=0.606833f, C2=0.960863f; break;
	  case 2: s1=0.0185603f, s2=0.0624481f, C1=0.553192f, C2=0.974129f; break;
	  }

	  xi = pf<C1 ? 1-s1*fabs(grnd(&s)) : pf<C2 ? 1+s2*log(xrnd(&s)) : -1.f;

	  if(xi<=-1.f) xi=2*sqrt(xrnd(&s))-1;
	  float si=sqrt(1-xi*xi); n=turn(xi, si, n, &s);
	}
	else if(p.ka!=0){ // old 2d gaussian
	  if(p.ka<0) xi=2*xrnd(&s)-1;
	  else do{ xi=1+p.ka*log(xrnd(&s)); } while (xi<-1);
	  float si=sqrt(1-xi*xi); n=turn(xi, si, n, &s);
	}

	if(p.ka>0){
	  float b=r.x*n.x+r.y*n.y+r.z*n.z;
	  float c=FLZ*FLZ+FLR*FLR-OMR*OMR;
	  float t=sqrt(b*b-c)-b;
	  r+=t*(float4)(n.xyz, e.ocv);
	  if(fabs(r.z)<FLB) ofla=-2;
	  else if(r.z<0) if(xrnd(&s)<0.9) ofla=-3;

	  if(p.ka>1050){
	    float rc=0.023f;
	    float dc=OMR+rc;

	    float ph=radians(p.ka-1080);
	    {
	      float drx=r.x-dc*cos(ph);
	      float dry=r.y-dc*sin(ph);
	      float a=n.x*n.x+n.y*n.y;
	      if(a>0){
		float b=n.x*drx+n.y*dry;
		float c=drx*drx+dry*dry-rc*rc;
		float D=b*b-a*c;
		if(D>=0){
		  float h1=(-b+sqrt(D))/a;
		  if(h1>0) ofla=-4;
		}
	      }
	    }
	  }

	  r+=p.r;
	}
      }
      else{
	float xi;
	if(p.ka>0){
	  if(p.type==5){
	    xi=2*sqrt(xrnd(&s))-1;
	  }
	  else{
	    do{ xi=1+p.ka*log(xrnd(&s)); } while (xi<-1);
	  }
	  float si=sqrt(1-xi*xi); n=turn(xi, si, n, &s);
	}
      }
    }
    else{
      fla=-1, ofla=-1;
      if(l>0) l*=xrnd(&s);
      else l=p.b*mrnd(p.a, &s);

      if(l>0){
	r.xyz+=l*n.xyz;
	r.w+=l*e.ocv/p.beta;
      }

      if(p.tau>0){ // isotropic delayed emmission
	r.w-=p.tau*log(xrnd(&s));
	float cs=xrnd(&s);
	float si=sqrt(1-cs*cs);
	n=turn(cs, si, n, &s);
      }
      else{
	float cs=w->coschr, si=w->sinchr;

	if(p.f<xrnd(&s)){ // cascade particle directions
	  const float a=0.39f, b=2.61f;
	  const float I=1-exp(-b*exp2(a));
	  float cs=max(1-pow(-log(1-xrnd(&s)*I)/b, 1/a), -1.0f);
	  float si=sqrt(1-cs*cs); n=turn(cs, si, n, &s);
	}
	else{
	  float beta=p.beta;
	  if(p.type<0){ // muon end point; assuming p.beta=1
	    const float ar=0.26f*0.9216f/0.105658389f;  // a [GeV/mwe] * density [mwe/m] / muon rest mass [GeV]
	    float dx=ar*(p.n.w-l);
	    beta=sqrt(dx*(2+dx))/(1+dx);
	    if(beta>=cs) r.w-=e.ocv*(sqrt(dx*dx+1)-asinh(1/dx)-dx)/ar;
	  }
	  if(beta<1){
	    float xi=cs/beta;
	    if(xi<=1){
	      float sx=sqrt(1-xi*xi);
	      if(p.type<0) if(sx<xrnd(&s)*si) ofla=-3;
	      cs=xi, si=sx;
	    }
	    else ofla=-2;
	  }
	}
	n=turn(cs, si, n, &s); // sampling cherenkov cone
      }
    }

    pbuf f; f.r=r, f.n=n; f.q=j; f.i=niw; f.fla=fla, f.ofla=ofla; bf[i]=f;
  }
  barrier(CLK_GLOBAL_MEM_FENCE);

  int ofla=-1;
  for(uint i=idx; i<num; TOT==0 && (i=atomic_add(&e.hidx, get_num_groups(0)))){
    int om=-1;
    if(TOT==0){ // initialize photon
      pbuf f=bf[i];
      r=f.r; n=f.n;
      w=&ez->w[f.q];
      niw=f.i;
      om=f.fla;
      ofla=f.ofla;

      TOT=-log(xrnd(&s)), SCA=0; // TOT = random number sampled from exponential to give distance to absorption
      if(ofla<-1) TOT=0;
    }
    float sqd;

    if(SCA==0) SCA=-log(xrnd(&s)), old=om; // SCA = random number sampled from exponential to give distance to next scattering
    float sca, tot; // get distance for overburden

    int sign=n.z<0?-1:1;
    float anz=sign*n.z;

    float edh = e.dh;
    float z = r.z - zshift(&e, r, e.vthk==0?0:sign, &edh, ez); // tilt correction

    float nr=1.f;
    int I, J;

    {
      z=(z-e.hmin)*e.rdh;
      I=clamp(convert_int_sat_rte(z), 0, e.size-1); // curent layer

      if(TOT>0){ // anisotropic absorption
	float n1= e.azx*n.x+e.azy*n.y;
	float n2=-e.azy*n.x+e.azx*n.y;
	float n3= n.z;

	if(e.fr!=0){
	  aniz az=e.az[I];
	  float k1=az.k1, k2=az.k2;
	  float kz=1/(k1*k2);

	  float s1=n1*n1, l1=k1*k1;
	  float s2=n2*n2, l2=k2*k2;
	  float s3=n3*n3, l3=kz*kz;

	  float B2=nr/l1+nr/l2+nr/l3;
	  float nB=s1/l1+s2/l2+s3/l3;
	  float An=s1*l1+s2*l2+s3*l3;

	  nr=pow((B2-nB)*An/2, e.fr);
	}

	{ // new absorption anisotropy
	  n1*=e.k1; n2*=e.k2; n3*=e.kz;
	  nr*=rsqrt(n1*n1+n2*n2+n3*n3);
	}

	TOT/=nr;
      }

      float ahx=edh*(0.5f-sign*(z-I)); // next layer boundary depending on direction

      J=I; // step trough layers until next scattering / absorption
      sca=0, tot=0;
      float ais=SCA, aia=TOT;
      while(true){
	float dsca=ais/w->z[J].sca, dtot=aia/w->z[J].abs;
	float dist=min(dsca, dtot)*anz;
	if((sign<0?J>0:J<e.size-1) && dist>ahx){
	  float del=ahx/anz;
	  sca+=del, tot+=del;
	  ais-=w->z[J].sca*del, aia-=w->z[J].abs*del;
	  ahx=edh;
	  J+=sign;
	}
	else{
	  sca+=dsca, tot+=dtot;
	  break;
	}
      }

      if(sca<0) sca=0;
      if(tot<0) tot=0;

      // get overburden for distance
      if(tot<sca) sca=tot, tot=0; else tot=(tot-sca)*w->z[J].abs;
   }

    om=-1;
    float del=sca;
    float hi=sca, hf=0;
    { // sphere
      float2 ri = r.xy, rf = r.xy + sca*n.xy;
      float2 pi = ctr(&e, ri), pf = ctr(&e, rf);
      ri = min(pi, pf)-e.rx, rf = max(pi, pf)+e.rx;

      int2 xl = (int2)(clamp(convert_int_sat_rte((ri.x-e.cl[0])*e.crst[0]), 0, e.cn[0]),
		       clamp(convert_int_sat_rte((ri.y-e.cl[1])*e.crst[1]), 0, e.cn[1]));

      int2 xh = (int2)(clamp(convert_int_sat_rte((rf.x-e.cl[0])*e.crst[0]), -1, e.cn[0]-1),
		       clamp(convert_int_sat_rte((rf.y-e.cl[1])*e.crst[1]), -1, e.cn[1]-1));

      for(int i=xl.x, j=xl.y; i<=xh.x && j<=xh.y; ++j<=xh.y?:(j=xl.y,i++)) for(ushort k=e.is[i][j]; k!=0x8000; ){
	ushort m=e.ls[k];
	__local line * s = & e.sc[m&0x7fff];
	k=m&0x8000?0x8000:k+1;

	float b=0, c=0, dr;
	dr=s->x-r.x;
	b+=n.x*dr; c+=dr*dr;
	dr=s->y-r.y;
	b+=n.y*dr; c+=dr*dr;

	float np=1-n.z*n.z;
	float D=b*b-(c-s->r*s->r)*np;
	if(D>=0){
	  D=sqrt(D);
	  float h1=b-D, h2=b+D;
	  if(h2>=0 && h1<=sca*np){
	    if(np>XXX){
	      h1/=np, h2/=np;
	      if(h1<0) h1=0; if(h2>sca) h2=sca;
	    }
	    else h1=0, h2=sca;
	    h1=r.z+n.z*h1, h2=r.z+n.z*h2;
	    float zl, zh;
	    if(n.z>0) zl=h1, zh=h2;
	    else zl=h2, zh=h1;

	    int omin=0, omax=s->max;
	    int n1=s->n-omin+clamp(convert_int_sat_rtp(omin-(zh-s->dl-s->h)*s->d), omin, omax+1);
	    int n2=s->n-omin+clamp(convert_int_sat_rtn(omin-(zl-s->dh-s->h)*s->d), omin-1, omax);

	    for(int l=n1; l<=n2; l++) if(l!=old){
	      // if(l==ofla) continue;
	      DOM dom = oms[l];

	      float a=0, b=0, c=0, dr;
	      float f=dom.F>0?square(1/dom.F):0;

	      dr=dom.r[0]-r.x; a+=  n.x*n.x; b+=  n.x*dr; c+=  dr*dr;
	      dr=dom.r[1]-r.y; a+=  n.y*n.y; b+=  n.y*dr; c+=  dr*dr;
	      dr=dom.r[2]-r.z; a+=f*n.z*n.z; b+=f*n.z*dr; c+=f*dr*dr;

	      if(a>0){
		b/=a;

		float xR=(l==ofla?1:e.xR);
		float D=b*b-(c-square(xR*dom.R))/a;
		if(D>=0){
		  sqd=sqrt(D);
		  float h=b-sqd;

		  if(dom.F>0 || fabs(dr-n.z*h)<-xR*dom.F*dom.R){
		    sqd/=xR; h=b-sqd;
		    if(h>0 && h<=del) om=l, del=h+0.0f;
		  }
		}
	      }
	    }
	  }
	  if(e.hr>0 && (e.hifl==0 || (e.hifl<0 && r.w<-e.hifl) || (e.hifl>0 && r.w>e.hifl))){
	    float D=b*b-(c-e.hr2)*np;
	    if(D>0){
	      D=sqrt(D);
	      float h1=b-D, h2=b+D;
	      if(h2>=0 && h1<=sca*np){
		if(np>XXX){
		  h1/=np, h2/=np;
		  if(h1<0) h1=0; if(h2>sca) h2=sca;
		}
		else h1=0, h2=sca;
		if(h1<hi && h2>sqrt(XXX)*e.hr) hi=h1, hf=h2;
	      }
	    }
	  }
	}
      }
    }

    float fin=min(del, hi);
    bool hole=fin<sca;

    if(hole){ // hole ice propagation code
      { // get overburden for distance
	float xs=0, xa=0;

	float y=z+n.z*fin*e.rdh;
	J=clamp(convert_int_sat_rte(y), 0, e.size-1);

	if(I==J) xs=fin*w->z[I].sca, xa=fin*w->z[I].abs;
	else{
	  float h=0.5f-sign*(z-I);
	  float g=0.5f+sign*(y-J);
	  xs=h*w->z[I].sca+g*w->z[J].sca;
	  xa=h*w->z[I].abs+g*w->z[J].abs;
	  for(int i=I+sign; i!=J; i+=sign){
	    xs+=w->z[i].sca;
	    xa+=w->z[i].abs;
	  }
	  float f=edh/anz;
	  xs*=f, xa*=f;
	}
	if(xs<0) xs=0;
	if(xa<0) xa=0;

	SCA-=xs, TOT-=xa;
      }
      TOT*=nr;

      if(hi<del){
	float fin=min(hi+min(SCA/e.hs, TOT/e.ha), hf);

	if(om!=-1){
	  if(fin<del) del=fin, om=-1;
	}
	else{
	  if(fin<del) del=fin;
	  else del=min(fin, sca);
	}

	fin=del-hi; SCA-=fin*e.hs, TOT-=fin*e.ha;
      }
    }
    else SCA=0, TOT=tot*nr;

    { // advance
      r+=del*n;
    }

    if(om!=-1){ // DOM collision was detected
      if(om==ofla) TOT=0;
      else{
	hit h; h.i=om; h.t=r.w; h.n=niw; h.z=w->wvl;

	h.pth=acos(n.z); h.pph=atan2(n.y, n.x);

	sqd*=e.xR-1;
	const DOM dom = oms[om];
	float dx=dom.r[0]-r.x+sqd*n.x;
	float dy=dom.r[1]-r.y+sqd*n.y;
	float dz=dom.r[2]-r.z+sqd*n.z;

	h.dth=acos(clamp(dz/(e.xR*dom.R*dom.F), -1.f, 1.f)); h.dph=atan2(dy, dx);

	{
	  uint j = atomic_inc(&ed->hidx);
	  if(j<e.hnum) eh[j]=h;
	}

	if(e.xR==1) TOT=0; else old=om;
      }
    }

    if(TOT<XXX) TOT=0;
    else{
      if(e.sum>0 && fin>0){ // birefringence
	float ra; // e.az[J].ra*fin
	float rb; // e.az[J].rb*fin

	{ // get overburden for distance
	  if(I==J) ra=fin*e.az[I].ra, rb=fin*e.az[I].rb;
	  else{
	    float y=z+n.z*fin*e.rdh;
	    float h=0.5f-sign*(z-I);
	    float g=0.5f+sign*(y-J);
	    ra=h*e.az[I].ra+g*e.az[J].ra;
	    rb=h*e.az[I].rb+g*e.az[J].rb;
	    for(int i=I+sign; i!=J; i+=sign){
	      ra+=e.az[i].ra;
	      rb+=e.az[i].rb;
	    }
	    float f=edh/anz;
	    ra*=f, rb*=f;
	  }
	  if(ra<0) ra=0;
	  if(rb<0) rb=0;
	}

	float dot=e.azx*n.x+e.azy*n.y;
	float sdt=sqrt(1-dot*dot);
	float cdt=fabs(dot);

	float sx=sqrt(ra)*max(0.f, e.bfr[0]*exp(-e.bfr[1]*pow(atan(e.bfr[3]*sdt), e.bfr[2])));
	float sy=sqrt(ra)*max(0.f, e.bfr[4]*exp(-e.bfr[5]*pow(atan(e.bfr[7]*sdt), e.bfr[6])));
	float mx=rb*max(0.f, e.bfr[8]*atan(e.bfr[11]*sdt*cdt)*exp(-e.bfr[9]*sdt+e.bfr[10]*cdt));

	float3 flow=(float3)(e.azx, e.azy, 0);
	float3 qz=n.xyz;
	float3 qx=flow-dot*qz; if(dot<0) qx=-qx;
	float qn=qx.x*qx.x+qx.y*qx.y+qx.z*qx.z;
	if(qn>0) qx*=rsqrt(qn); else qx=(float3)(0, 0, 1);
	float3 qy=(float3)(qz.y*qx.z-qz.z*qx.y, qz.z*qx.x-qz.x*qx.z, qz.x*qx.y-qz.y*qx.x);

	float2 rnd=grnd2(&s);
	float dnx=sx*rnd.x+mx;
	float dny=sy*rnd.y;

	if(del<=hi){ // curve
	  rnd=grnd2(&s);
	  float vnx=dnx+sx*rnd.x*rsqrt(3.f);
	  float vny=dny+sy*rnd.y*rsqrt(3.f);

	  vnx/=2, vny/=2;
	  float den=vnx*vnx+vny*vny;
	  float dnz=-den; den=2/(1+den);
	  vnx*=den, vny*=den, dnz*=den;
	  float3 m=n.xyz+qx*vnx+qy*vny+qz*dnz; m=mx_normalize(m);

	  del-=fin/2;
	  r.xyz+=del*(m-n.xyz);
	}

	dnx/=2, dny/=2;
	float den=dnx*dnx+dny*dny;
	float dnz=-den; den=2/(1+den);
	dnx*=den, dny*=den, dnz*=den;
	n.xyz+=qx*dnx+qy*dny+qz*dnz; n=my_normalize(n);
      }

      if(SCA<XXX){
	SCA=0;
	float sf, g, gr;
	if(hole){
	  sf=e.SF, g=e.G, gr=e.GR;
	}
	else{
	  sf=e.sf, g=e.g, gr=e.gr;
	}

	float xi=xrnd(&s);
	if(xi>sf){
	  xi=(1-xi)/(1-sf);
	  xi=2*xi-1;
	  if(g!=0){
	    float ga=(1-g*g)/(1+g*xi);
	    xi=(1+g*g-ga*ga)/(2*g);
	  }
	}
	else{
	  xi/=sf;
	  xi=2*pow(xi, gr)-1;
	}

	if(xi>1) xi=1; else if(xi<-1) xi=-1;

	aniz az=e.az[J];
	float k1=az.k1, k2=az.k2;
	float okz=k1*k2;

	if(!hole){ // if not in hole ice, rotate coordinate system for anisotropic scattering
	  float n1=( e.azx*n.x+e.azy*n.y)*k1;
	  float n2=(-e.azy*n.x+e.azx*n.y)*k2;
	  n.x=n1*e.azx-n2*e.azy;
	  n.y=n1*e.azy+n2*e.azx;
	  n.z/=okz;
	  n=my_normalize(n);
	}

	float si=sqrt(1-xi*xi); // perform scattering
	n=turn(xi, si, n, &s);

	if(!hole){ // rotate back into default coordinate system
	  float n1=( e.azx*n.x+e.azy*n.y)/k1;
	  float n2=(-e.azy*n.x+e.azx*n.y)/k2;
	  n.x=n1*e.azx-n2*e.azy;
	  n.y=n1*e.azy+n2*e.azx;
	  n.z*=okz;
	  n=my_normalize(n);
	}
      }

      if(!isfinite(TOT+SCA)) TOT=0, atomic_inc(&ed->ab);
    }
  }

  {
    ez->rs[s.w]=s.x | (ulong) s.y << 32;
    barrier(CLK_LOCAL_MEM_FENCE);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
  }

}

#undef sin
#undef cos
#undef pow
#undef exp
#undef log
#undef exp2
#undef sqrt
#undef rsqrt

));
#endif
