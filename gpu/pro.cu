#ifdef XCPU
#define __device__
#define __global__

#define rsqrtf 1/sqrtf
#define __float2int_rn (int)lroundf
#define __float2int_ru (int)ceilf
#define __float2int_rd (int)floorf

struct int2{
  int x, y;
};

struct int3:int2{
  int z;
};

struct uint4{
  unsigned int x, y, z, w;
};

float __int_as_float(unsigned int x){
  union{
    unsigned int i;
    float f;
  };
  i=x; return f;
}

unsigned int atomicAdd(unsigned int * i, unsigned int j){
  unsigned int k=*i; *i+=j;
  return k;
}

struct ThreadIdx{
  int x;
} threadIdx;

struct BlockDim{
  int x;
} blockDim;

unsigned int seed=0;

#if defined(__APPLE_CC__) || defined(__FreeBSD__)
void sincosf(float x, float * s, float * c){ *s = sin(x); *c = cos(x); }
#endif
#endif

__device__ float xrnd(uint4 & s){
  unsigned long long sda = s.y;
#ifdef XCPU
  sda += s.z * (unsigned long long) s.x;
#else
  asm("mad.wide.u32 %0, %1, %2, %0;" : "+l"(sda) : "r"(s.x), "r"(s.z));
#endif
  s.x = sda; s.y = sda >> 32;
  unsigned int tmp = s.x >> 9;
  return 2.0f-__int_as_float(tmp|0x3f800000);
}

__device__ float mrnd(float k, uint4 & s){  // gamma distribution
  float x;
  if(k<1){  // Weibull algorithm
    float c=1/k;
    float d=(1-k)*powf(k, 1/(c-1));
    float z, e;
    do{
      z=-logf(xrnd(s));
      e=-logf(xrnd(s));
      x=powf(z, c);
    } while(z+e<d+x);
  }
  else{  // Cheng's algorithm
    float b=k-logf(4.0f);
    float l=sqrtf(2*k-1);
    float c=1+logf(4.5f);
    float u, v, y, z, r;
    do{
      u=xrnd(s); v=xrnd(s);
      y=-logf(1/v-1)/l;
      x=k*expf(y);
      z=u*v*v;
      r=b+(k+l)*y-x;
    } while(r<4.5f*z-c && r<logf(z));
  }
  return x;
}

__device__ float grnd(uint4 & s){  // gaussian distribution
  return sqrtf(-2*logf(xrnd(s)))*sinf(2*FPI*xrnd(s));
}

__device__ void grnd2(float2 & p, uint4 & s){  // gaussian distribution
  float r=sqrtf(-2*logf(xrnd(s)));
  sincosf(2*FPI*xrnd(s), &p.y, &p.x);
  p.x*=r, p.y*=r;
}

__device__ void mx_normalize(float3 & n){
  float r=rsqrtf(n.x*n.x+n.y*n.y+n.z*n.z);
  n.x*=r, n.y*=r, n.z*=r;
}

__device__ void my_normalize(float4 & n){
  float r=rsqrtf(n.x*n.x+n.y*n.y+n.z*n.z);
  n.x*=r, n.y*=r, n.z*=r;
}

__device__ void swap(float & x, float & y){
  float a=x; x=y; y=a;
}

__device__ void rotate(float & cs, float & si, float4 & n, uint4 & s){
  float3 p1, p2;
  int i=0;
  {
    float3 r;
    r.x=n.x*n.x, r.y=n.y*n.y, r.z=n.z*n.z;
    if(r.y>r.z){
      if(r.y>r.x) i=(swap(n.x,n.y),swap(r.x,r.y),1);
    }
    else{
      if(r.z>r.x) i=(swap(n.x,n.z),swap(r.x,r.z),2);
    }

    r.y=rsqrtf(r.x+r.y); p1.x=-n.y*r.y; p1.y=n.x*r.y; p1.z=0;
    r.z=rsqrtf(r.x+r.z); p2.x=-n.z*r.z; p2.y=0; p2.z=n.x*r.z;
  }

  {
    float4 q1;

    q1.x=p1.x-p2.x; q1.y=p1.y-p2.y; q1.z=p1.z-p2.z;
    p2.x+=p1.x; p2.y+=p1.y; p2.z+=p1.z;

    q1.w=rsqrtf(q1.x*q1.x+q1.y*q1.y+q1.z*q1.z);
    p1.x=q1.x*q1.w; p1.y=q1.y*q1.w; p1.z=q1.z*q1.w;

    q1.w=rsqrtf(p2.x*p2.x+p2.y*p2.y+p2.z*p2.z);
    p2.x*=q1.w; p2.y*=q1.w; p2.z*=q1.w;
  }

  {
    float2 p;
    float xi=2*FPI*xrnd(s);
    sincosf(xi, &p.y, &p.x);

    n.x=cs*n.x+si*(p.x*p1.x+p.y*p2.x);
    n.y=cs*n.y+si*(p.x*p1.y+p.y*p2.y);
    n.z=cs*n.z+si*(p.x*p1.z+p.y*p2.z);

    my_normalize(n);
    if(i==1) swap(n.x,n.y); else if(i==2) swap(n.x,n.z);
  }
}

#ifndef XCPU
__host__ __device__ int __float2int_rd(float x);
__host__ int __float2int_rd(float x){ return (int)floorf(x); }
__host__
#endif

__device__ float square(float x){
  return x*x;
}

#ifndef XCPU
__host__
#endif

__device__ int findcol(dats & d, int x, int y){
  return x>=0 && x<d.mnum[0] && y>=0 && y<d.mnum[1] ? d.mcol[x][y] : -1;
}

#ifndef XCPU
__host__
#endif

__device__ float midp(float2 q, float z){
  return q.x+q.y*z;
}

#ifndef XCPU
__host__
#endif

__device__ float interpolate(dats & d, int j, float nr, float z){
  int i=j-1, k=min(max(__float2int_rd(z), 0), d.lpts);
  datz & p = * d.z;
  return ( midp(p.lp[j][k], z-k)*(nr-d.lr[i]) +
	   midp(p.lp[i][k], z-k)*(d.lr[j]-nr) )/(d.lr[j]-d.lr[i]);
}

#ifndef XCPU
__host__
#endif

__device__ float interpolate(dats & d, int3 & c, float3 & w, float z){
  int k=min(max(__float2int_rd(z), 0), d.lpts);
  datz & p = * d.z;
  return
    ( c.x<0?0: midp(p.lp[c.x][k], z-k)*w.x )+
    ( c.y<0?0: midp(p.lp[c.y][k], z-k)*w.y )+
    ( c.z<0?0: midp(p.lp[c.z][k], z-k)*w.z );
}

#ifndef XCPU
__host__
#endif

__device__ float zshift(dats & d, float4 & r, int sign, float & dz){
  if(d.tmod==1){
    float nr=d.lnx*r.x+d.lny*r.y-d.r0;
    int j=1; for(; j<LMAX; j++) if(nr<d.lr[j] || j==d.lnum-1) break;

    float z1=(r.z-d.lmin)*d.lrdz;
    float r1=interpolate(d, j, nr, z1);
    if(sign!=0){
      float z2=z1+sign*dz*d.lrdz;
      float r2=interpolate(d, j, nr, z2);
      dz /= 1-sign*(r2-r1)/ dz;
    }
    return r1;
  }
  else if(d.tmod==2){
    float2 q; // ctr(d, r, q);
    q.x=d.cb[0][0]*r.x+d.cb[1][0]*r.y;
    q.y=d.cb[0][1]*r.x+d.cb[1][1]*r.y;

    q.x=(q.x-d.mmin[0])/d.mstp[0];
    q.y=(q.y-d.mmin[1])/d.mstp[1];

    { // extrapolate
      float x=q.x, y=q.y, xy;
      int n0=d.mnum[0]-1, n1=d.mnum[1]-1, c0=d.mcut[0]-1, c1=d.mcut[1]-1;

      if(y<0) if(x-y>0 && x<=c0){ x=(x-y)*c0/(c0-y); y=0; } // 1
      if(x<0) if(y-x>=0 && y<c1){ y=(y-x)*c1/(c1-x); x=0; } // 2

      if(y>n1) if(x>=n1-c1 && y-x>n1-n0){ x=n1-c1+(n0-(n1-c1))*(x-(n1-c1))/(y-(n1-n0)-(n1-c1)); y=n1; } // 5
      if(x>n0) if(y>n0-c0 && y-x<=n1-n0){ y=n0-c0+(n1-(n0-c0))*(y-(n0-c0))/(x-(n0-n1)-(n0-c0)); x=n0; } // 6

      if(y-x>c1) if(y>=c1 && x<n1-c1){ xy=y-x-c1; xy=c1+(x+y-c1+xy)*(n1-c1)/(n1-c1+xy); x=(xy-c1)/2; y=(xy+c1)/2; } // 3
      if(x-y>c0) if(x>c0 && y<=n0-c0){ xy=x-y-c0; xy=c0+(x+y-c0+xy)*(n0-c0)/(n0-c0+xy); x=(xy+c0)/2; y=(xy-c0)/2; } // 4

      q.x=x, q.y=y;
    }

    int2 qn;
    qn.x=min(max(__float2int_rd(q.x), 0), d.mnum[0]);
    qn.y=min(max(__float2int_rd(q.y), 0), d.mnum[1]);

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

    float z1=(r.z-d.lmin)*d.lrdz;
    float r1=interpolate(d, c, w, z1);
    if(sign!=0){
      float z2=z1+sign*dz*d.lrdz;
      float r2=interpolate(d, c, w, z2);
      dz /= 1-sign*(r2-r1)/ dz;
    }

    return r1;
  }
  else return 0;
}

__device__ void ctr(dats & d, float2 & r, float2 & p){
  p.x=d.cb[0][0]*r.x+d.cb[1][0]*r.y;
  p.y=d.cb[0][1]*r.x+d.cb[1][1]*r.y;
}

#ifndef XCPU
__device__ inline unsigned int smid(){
  unsigned int r;
  asm volatile("mov.u32 %0, %%smid;" : "=r"(r));
  return r;
}
#endif

#ifdef XCPU
#define XINC i+=eidx
#define XIDX e.gridDim*blockDim.x
#else
#define XINC i=atomicAdd(&eidx, e.gridDim)
#define XIDX e.gridDim*blockDim.x+e.blockIdx
#endif

__global__ void propagate(dats * ed, unsigned int num){
  uint4 s;
  unsigned int niw=0;
#ifdef XCPU
  float4 n;
  float4 r;
  dats & e = * ed;
  static unsigned int eidx;
  if(threadIdx.x==0) eidx = XIDX;
#else
  float4 n={0,0,0,0};
  float4 r={0,0,0,0};
  __shared__ dats e;
  unsigned int & eidx = e.hidx;

  if(num==0){
    ed->hidx=0;
    ed->tn=-1U;
    ed->tx=0;
    ed->ab=0;
    ed->mp=0;
    __threadfence();
    return;
  }

  if(threadIdx.x==0){
    e=*ed; e.tn=clock();
    e.blockIdx=smid()==e.blockIdx?-1:(int)atomicAdd(&ed->mp, 1);
    eidx=XIDX;
  }
  __syncthreads();

  if(e.blockIdx==-1) return;
#endif

  ices * w;
  const unsigned int idx=threadIdx.x*e.gridDim+e.blockIdx;

  {
#ifndef XCPU
    const unsigned int & seed = idx;
#endif
    s.w=seed%e.rsize;
    s.x=e.z->rs[s.w];
    s.y=e.z->rs[s.w] >> 32;
    s.z=e.z->rm[s.w];
  }

  int old;
  float TOT=0, SCA;

  for(unsigned int i=idx, pj=-1U, pn=0; i<num; i+=e.gridDim*blockDim.x){
    while(e.gini+i%e.gspc+(i/e.gspc)*e.gtot>=pn) pn+=e.pz[++pj].num;

    unsigned int j=min(__float2int_rd(WNUM*xrnd(s)), WNUM-1);
    w=&e.z->w[j];
    n.w=w->ocm;

    photon p=e.pz[pj];
    r=p.r, n.x=p.n.x, n.y=p.n.y, n.z=p.n.z;
    float l=p.n.w; niw=p.q;
    int fla, ofla;

    if(p.type>0){
      fla=p.fla, ofla=p.ofla;

      if(p.type<=4){
	float xi=xrnd(s);
	if(p.fldr<0) xi*=2*FPI;
	else{
	  int r=__float2int_rd(p.fldr/360)+1;
	  int s=__float2int_rd(xi*r);
	  xi=(p.fldr+s*360/r)*fcv;
	}
	sincosf(xi, &n.y, &n.x);

	float FLZ=0.07735f, FLR=0.119f-0.008f, FLB=0.03396f;
	if(p.ka>0) r.x=FLR*n.x, r.y=FLR*n.y, r.z=FLZ, r.w=0;

	sincosf(p.up, &n.z, &xi);
	n.x*=xi; n.y*=xi;

	if(p.ka>999.f){
	  float pf=xrnd(s);

	  float   s1=0.0191329f, s2=0.0686944f, C1=0.583583f, C2=0.967762f;
	  switch(p.type){
	  case 1: s1=0.0196242f, s2=0.0750278f, C1=0.606833f, C2=0.960863f; break;
	  case 2: s1=0.0185603f, s2=0.0624481f, C1=0.553192f, C2=0.974129f; break;
	  }

	  xi = pf<C1 ? 1-s1*fabsf(grnd(s)) : pf<C2 ? 1+s2*logf(xrnd(s)) : -1.f;

	  if(xi<=-1.f) xi=2*sqrtf(xrnd(s))-1;
	  float si=sqrtf(1-xi*xi); rotate(xi, si, n, s);
	}
	else if(p.ka!=0){ // old 2d gaussian
	  if(p.ka<0) xi=2*xrnd(s)-1;
	  else do{ xi=1+p.ka*logf(xrnd(s)); } while (xi<-1);
	  float si=sqrtf(1-xi*xi); rotate(xi, si, n, s);
	}

	if(p.ka>0){
	  float b=r.x*n.x+r.y*n.y+r.z*n.z;
	  float c=FLZ*FLZ+FLR*FLR-OMR*OMR;
	  float t=sqrtf(b*b-c)-b;
	  r.x+=t*n.x, r.y+=t*n.y, r.z+=t*n.z, r.w+=t*e.ocv;
	  if(fabsf(r.z)<FLB) ofla=-2;
	  else if(r.z<0) if(xrnd(s)<0.9) ofla=-3;

	  if(p.ka>1050){
	    float rc=0.023f;
	    float dc=OMR+rc;

	    float ph=(p.ka-1080)*fcv;
	    {
	      float drx=r.x-dc*cosf(ph);
	      float dry=r.y-dc*sinf(ph);
	      float a=n.x*n.x+n.y*n.y;
	      if(a>0){
		float b=n.x*drx+n.y*dry;
		float c=drx*drx+dry*dry-rc*rc;
		float D=b*b-a*c;
		if(D>=0){
		  float h1=(-b+sqrtf(D))/a;
		  if(h1>0) ofla=-4;
		}
	      }
	    }
	  }

	  r.x+=p.r.x, r.y+=p.r.y, r.z+=p.r.z, r.w+=p.r.w;
	}
      }
      else{
	float xi;
	if(p.ka>0){
	  if(p.type==5){
	    xi=2*sqrtf(xrnd(s))-1;
	  }
	  else{
	    do{ xi=1+p.ka*logf(xrnd(s)); } while (xi<-1);
	  }
	  float si=sqrtf(1-xi*xi); rotate(xi, si, n, s);
	}
      }
    }
    else{
      fla=-1, ofla=-1;
      if(l>0) l*=xrnd(s);
      else l=p.b*mrnd(p.a, s);

      if(l>0){
	r.x+=l*n.x, r.y+=l*n.y, r.z+=l*n.z;
	r.w+=l*e.ocv/p.beta;
      }

      if(p.tau>0){ // isotropic delayed emmission
	r.w-=p.tau*logf(xrnd(s));
	float cs=xrnd(s);
	float si=sqrtf(1-cs*cs);
	rotate(cs, si, n, s);
      }
      else{
	float cs=w->coschr, si=w->sinchr;

	if(p.f<xrnd(s)){ // cascade particle directions
	  const float a=0.39f, b=2.61f;
	  const float I=1-expf(-b*exp2(a));
	  float cs=max(1-powf(-logf(1-xrnd(s)*I)/b, 1/a), -1.0f);
	  float si=sqrtf(1-cs*cs); rotate(cs, si, n, s);
	}
	else{
	  float beta=p.beta;
	  if(p.type<0){ // muon end point; assuming p.beta=1
	    const float ar=0.26f*0.9216f/0.105658389f;  // a [GeV/mwe] * density [mwe/m] / muon rest mass [GeV]
	    float dx=ar*(p.n.w-l);
	    beta=sqrtf(dx*(2+dx))/(1+dx);
	    if(beta>=cs) r.w-=e.ocv*(sqrtf(dx*dx+1)-asinhf(1/dx)-dx)/ar;
	  }
	  if(beta<1){
	    float xi=cs/beta;
	    if(xi<=1){
	      float sx=sqrtf(1-xi*xi);
	      if(p.type<0) if(sx<xrnd(s)*si) ofla=-3;
	      cs=xi, si=sx;
	    }
	    else ofla=-2;
	  }
	}
	rotate(cs, si, n, s); // sampling cherenkov cone
      }
    }

    pbuf f; f.r=r, f.n=n; f.q=j; f.i=niw; f.fla=fla, f.ofla=ofla; e.bf[i]=f;
  }
#ifndef XCPU
  __syncthreads();
#endif

  int ofla=-1;
  for(unsigned int i=idx; i<num; TOT==0 && (XINC)){
    int om=-1;
    if(TOT==0){ // initialize photon
      pbuf f=e.bf[i];
      r=f.r; n=f.n;
      w=&e.z->w[f.q];
      niw=f.i;
      om=f.fla;
      ofla=f.ofla;

      TOT=-logf(xrnd(s)), SCA=0; // TOT = random number sampled from exponential to give distance to absorption
      if(ofla<-1) TOT=0;
    }
    float sqd;

    if(SCA==0) SCA=-logf(xrnd(s)), old=om; // SCA = random number sampled from exponential to give distance to next scattering
    float sca, tot; // get distance for overburden

    int sign=n.z<0?-1:1;
    float anz=sign*n.z;

    float edh = e.dh;
    float z = r.z - zshift(e, r, e.vthk==0?0:sign, edh); // tilt correction

    float nr=1.f;
    int I, J;

    {
      z=(z-e.hmin)*e.rdh;
      I=min(max(__float2int_rn(z), 0), e.size-1); // curent layer

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

	  nr=powf((B2-nB)*An/2, e.fr);
	}

	{ // new absorption anisotropy
	  n1*=e.k1; n2*=e.k2; n3*=e.kz;
	  nr*=rsqrtf(n1*n1+n2*n2+n3*n3);
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
      float2 ri, rf, pi, pf;

      ri.x=r.x; rf.x=r.x+sca*n.x;
      ri.y=r.y; rf.y=r.y+sca*n.y;

      ctr(e, ri, pi); ctr(e, rf, pf);

      ri.x=min(pi.x, pf.x)-e.rx; rf.x=max(pi.x, pf.x)+e.rx;
      ri.y=min(pi.y, pf.y)-e.rx; rf.y=max(pi.y, pf.y)+e.rx;

      int2 xl, xh;

      xl.x=min(max(__float2int_rn((ri.x-e.cl[0])*e.crst[0]), 0), e.cn[0]);
      xh.x=max(min(__float2int_rn((rf.x-e.cl[0])*e.crst[0]), e.cn[0]-1), -1);

      xl.y=min(max(__float2int_rn((ri.y-e.cl[1])*e.crst[1]), 0), e.cn[1]);
      xh.y=max(min(__float2int_rn((rf.y-e.cl[1])*e.crst[1]), e.cn[1]-1), -1);

      for(int i=xl.x, j=xl.y; i<=xh.x && j<=xh.y; ++j<=xh.y?:(j=xl.y,i++)) for(unsigned short k=e.is[i][j]; k!=0x8000; ){
	unsigned short m=e.ls[k];
	line & s = e.sc[m&0x7fff];
	k=m&0x8000?0x8000:k+1;

	float b=0, c=0, dr;
	dr=s.x-r.x;
	b+=n.x*dr; c+=dr*dr;
	dr=s.y-r.y;
	b+=n.y*dr; c+=dr*dr;

	float np=1-n.z*n.z;
	float D=b*b-(c-s.r*s.r)*np;
	if(D>=0){
	  D=sqrtf(D);
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

	    int omin=0, omax=s.max;
	    int n1=s.n-omin+min(omax+1, max(omin, __float2int_ru(omin-(zh-s.dl-s.h)*s.d)));
	    int n2=s.n-omin+max(omin-1, min(omax, __float2int_rd(omin-(zl-s.dh-s.h)*s.d)));

	    for(int l=n1; l<=n2; l++) if(l!=old){
	      // if(l==e.fla) continue;
	      const DOM & dom=e.oms[l];

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
		  sqd=sqrtf(D);
		  float h=b-sqd;

		  if(dom.F>0 || fabsf(dr-n.z*h)<-xR*dom.F*dom.R){
		    sqd/=xR; h=b-sqd;
		    if(h>0 && h<=del) om=l, del=h;
		  }
		}
	      }
	    }
	  }
	  if(e.hr>0 && (e.hifl==0 || (e.hifl<0 && r.w<-e.hifl) || (e.hifl>0 && r.w>e.hifl))){
	    float D=b*b-(c-e.hr2)*np;
	    if(D>0){
	      D=sqrtf(D);
	      float h1=b-D, h2=b+D;
	      if(h2>=0 && h1<=sca*np){
		if(np>XXX){
		  h1/=np, h2/=np;
		  if(h1<0) h1=0; if(h2>sca) h2=sca;
		}
		else h1=0, h2=sca;
		if(h1<hi && h2>sqrtf(XXX)*e.hr) hi=h1, hf=h2;
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
	J=min(max(__float2int_rn(y), 0), e.size-1);

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
      r.x+=del*n.x;
      r.y+=del*n.y;
      r.z+=del*n.z;
      r.w+=del*n.w;
    }

#ifndef XCPU
    if(!isfinite(TOT) || !isfinite(SCA)) ed->bmp[atomicAdd(&ed->ab, 1)%4]=smid(), TOT=0, om=-1;
#endif

    if(om!=-1){ // DOM collision was detected
      if(om==ofla) TOT=0;
      else{
	hit h; h.i=om; h.t=r.w; h.n=niw; h.z=w->wvl;

	h.pth=acosf(n.z); h.pph=atan2f(n.y, n.x);

	sqd*=e.xR-1;
	const DOM & dom=e.oms[om];
	float dx=dom.r[0]-r.x+sqd*n.x;
	float dy=dom.r[1]-r.y+sqd*n.y;
	float dz=dom.r[2]-r.z+sqd*n.z;

	h.dth=acosf(min(max(dz/(e.xR*dom.R*dom.F), -1.f), 1.f)); h.dph=atan2f(dy, dx);

	{
	  unsigned int j = atomicAdd(&ed->hidx, 1);
	  if(j<e.hnum) e.hits[j]=h;
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
	float sdt=sqrtf(1-dot*dot);
	float cdt=fabsf(dot);

	float sx=sqrtf(ra)*max(0.f, e.bfr[0]*expf(-e.bfr[1]*powf(atanf(e.bfr[3]*sdt), e.bfr[2])));
	float sy=sqrtf(ra)*max(0.f, e.bfr[4]*expf(-e.bfr[5]*powf(atanf(e.bfr[7]*sdt), e.bfr[6])));
	float mx=rb*max(0.f, e.bfr[8]*atanf(e.bfr[11]*sdt*cdt)*expf(-e.bfr[9]*sdt+e.bfr[10]*cdt));

	float qx=e.azx-dot*n.x, qy=e.azy-dot*n.y, qz=-dot*n.z;
	if(dot<0) qx=-qx, qy=-qy, qz=-qz;

	float qn=qx*qx+qy*qy+qz*qz;
	if(qn>0){
	  float f=rsqrtf(qn);
	  qx*=f, qy*=f, qz*=f;
	}
	else{
	  qx=0, qy=0, qz=1;
	}
	float px=n.y*qz-n.z*qy, py=n.z*qx-n.x*qz, pz=n.x*qy-n.y*qx;

	float2 rnd; grnd2(rnd, s);
	float dnx=sx*rnd.x+mx;
	float dny=sy*rnd.y;

	if(del<=hi){ // curve
	  grnd2(rnd, s);
	  float vnx=dnx+sx*rnd.x*rsqrtf(3.f);
	  float vny=dny+sy*rnd.y*rsqrtf(3.f);

	  vnx/=2, vny/=2;
	  float den=vnx*vnx+vny*vny;
	  float dnz=-den; den=2/(1+den);
	  vnx*=den, vny*=den, dnz*=den;

	  float3 m;
	  m.x=n.x+qx*vnx+px*vny+n.x*dnz;
	  m.y=n.y+qy*vnx+py*vny+n.y*dnz;
	  m.z=n.z+qz*vnx+pz*vny+n.z*dnz;
	  mx_normalize(m);

	  del-=fin/2;
	  r.x+=del*(m.x-n.x);
	  r.y+=del*(m.y-n.y);
	  r.z+=del*(m.z-n.z);
	}

	dnx/=2, dny/=2;
	float den=dnx*dnx+dny*dny;
	float dnz=-den; den=2/(1+den);
	dnx*=den, dny*=den, dnz*=den;

	n.x+=qx*dnx+px*dny+n.x*dnz;
	n.y+=qy*dnx+py*dny+n.y*dnz;
	n.z+=qz*dnx+pz*dny+n.z*dnz;
	my_normalize(n);
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

	float xi=xrnd(s);
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
	  xi=2*powf(xi, gr)-1;
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
	  my_normalize(n);
	}

	float si=sqrtf(1-xi*xi); // perform scattering
	rotate(xi, si, n, s);

	if(!hole){ // rotate back into default coordinate system
	  float n1=( e.azx*n.x+e.azy*n.y)/k1;
	  float n2=(-e.azy*n.x+e.azx*n.y)/k2;
	  n.x=n1*e.azx-n2*e.azy;
	  n.y=n1*e.azy+n2*e.azx;
	  n.z*=okz;
	  my_normalize(n);
	}
      }
    }
  }

  {
    e.z->rs[s.w]=s.x | (unsigned long long) s.y << 32;
#ifndef XCPU
    __syncthreads();
    if(threadIdx.x==0){
      e.tx=clock();
      atomicMin(&ed->tn, e.tx-e.tn);
      atomicMax(&ed->tx, e.tx-e.tn);
    }
    __threadfence();
#endif
  }

}
