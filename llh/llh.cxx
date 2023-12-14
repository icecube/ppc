#include <set>
#include <map>
#include <vector>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

using namespace std;

#include "ppc.h"
#include "nnls.cxx"
#include "lssl.cxx"

double fdur=70.14, tau=4.00;
bool norm=false;
bool fail=false;
bool fast=false;
bool oldf=true;
bool mlpd=false;
bool flasher=false;
int munf=0;
int srep=1, drep=1;

int cnum=1;
double cbin=2./cnum;

const int num=200;
const bool cutz=true;  // cut end bins
const double bin=25;
const double qmin=0.1;
double qsat=500; // 2500

bool sumw=false;
const double nsrt=500e-9;
vector<double> noise, noise_q;

int fsep=1, loop=1000;
double lbin=bin, nsq=nsrt*lbin*drep;
double rr=0, dR=10, dA=30, dI=0, LX=0;
const double cv=M_PI/180, sol=0.299792458;

const double OMR=0.16510; // DOM radius [m]

void step(){
  cerr<<"rr="<<rr<<" dR="<<dR<<" dA="<<dA<<" dI="<<dI<<endl;
}

double square(double x){
  return x*x;
}

template <int n> class V{
  double x[n];

public:
  V(){
    for(int i=0; i<n; i++) x[i]=0;
  }

  double & operator[](int i){
    return x[i];
  }
};

struct xrnd{
  gsl_rng * r;

  xrnd(){
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);
  }

  ~xrnd(){
    gsl_rng_free(r);
  }

  double u(){
    return gsl_rng_uniform(r);
  }
} xr;

double xrnd(){
  // return drand48();
  return xr.u();
}

double cab(double *, double *);

void cab_fdf(const gsl_vector *x, void *par, double *f, gsl_vector *df){
  const int m=2;
  double X[m], D[m];
  for(int i=0; i<m; i++) X[i]=gsl_vector_get(x, i);
  double X0=X[0];
  X[0]*=X0; // ensure t=X[0]>0
  double result=cab(X, D);

  if(f!=NULL) *f = result;
  if(df!=NULL){
    D[0]*=2*X0;
    for(int i=0; i<m; i++) gsl_vector_set(df, i, D[i]);
  }
}

void cab_df(const gsl_vector *x, void *par, gsl_vector *df){
  cab_fdf(x, par, NULL, df);
}

double cab_f(const gsl_vector *x, void *par){
  double result;
  cab_fdf(x, par, &result, NULL);
  return result;
}

struct cable{
  double nr, cs, sn, rn, rc, dc;
  double nx, ny, nz, rx, ry, rz;
  double xx;
  bool cylr;

  cable(){
    rc=0.023f;
    dc=OMR+rc;
    xx=rc/100;

    {
      char * bmp=getenv("CYLR"); cylr=false;
      if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)) cylr=true;
      cerr<<"Cable: "<<(cylr?"cylindrical":"curved")<<"."<<endl;
    }
  }

  bool cyl(const float h[], const float r[], float ph){
    { // cylinder
      float drx=r[0]-dc*cos(cv*ph);
      float dry=r[1]-dc*sin(cv*ph);
      float a=h[0]*h[0]+h[1]*h[1];
      if(a>0){
	float b=-2*(h[0]*drx+h[1]*dry);
	float c=drx*drx+dry*dry-rc*rc;
	float D=b*b-4*a*c;
	if(D>=0){
	  float h1=(-b+sqrt(D))/(2*a);
	  if(h1>0) return true;
	}
      }
    }
    return false;
  }

  bool cab(const float h[], const float r[], float ph){
    if(cylr) return cyl(h, r, ph);

    nx=h[0], ny=h[1], nz=h[2];
    rx=r[0], ry=r[1], rz=r[2];

    cs=dc*cos(cv*ph), sn=dc*sin(cv*ph);
    nr=nx*rx+ny*ry+nz*rz; rn=nx*cs+ny*sn;

    {
      const gsl_multimin_fdfminimizer_type *T;
      gsl_multimin_fdfminimizer *s;

      gsl_vector *x;
      gsl_multimin_function_fdf G;

      const int m=2;
      G.n = m;
      G.f = cab_f;
      G.df = cab_df;
      G.fdf = cab_fdf;
      G.params = NULL;

      x = gsl_vector_alloc(m);
      double t0=0, aux=1-nz*nz;
      if(aux>0) t0=(nr-rn-nz*rz)/aux; if(t0<xx) t0=xx;
      double xi[m]={sqrt(t0), rz-t0*nz};
      for(int i=0; i<m; i++) gsl_vector_set(x, i, xi[i]);

      switch(2){
      case 0: T = gsl_multimin_fdfminimizer_conjugate_fr; break;
      case 1: T = gsl_multimin_fdfminimizer_conjugate_pr; break;
      case 2: T = gsl_multimin_fdfminimizer_vector_bfgs2; break;
      default: T = gsl_multimin_fdfminimizer_steepest_descent;
      }

      s = gsl_multimin_fdfminimizer_alloc(T, m);

      gsl_multimin_fdfminimizer_set(s, &G, x, 0.01, 0.1);

      for(int iter=0; iter<25; iter++) if(gsl_multimin_fdfminimizer_iterate(s)!=GSL_SUCCESS) break;
      double F=sqrt(s->f);

      if(false) if(F<rc){
	double t=gsl_vector_get(s->x, 0); t*=t;
	cout<<(rx-nx*t)<<" "<<(ry-ny*t)<<" "<<(rz-nz*t)<<endl;
      }

      gsl_multimin_fdfminimizer_free(s);
      gsl_vector_free(x);
      if(F<rc) return true;
    }

    return false;
  }

  double cab(double *X, double *D){
    double t=X[0], xi=X[1];
    const double q=0.6f;

    double aux=xi/q;
    double fx=exp(-aux*aux);
    double gx=-2*fx*aux/q;
    D[0]=-0.5*(nr-rn*fx-nz*xi-t);
    double a1=rx-nx*t-cs*fx;
    double a2=ry-ny*t-sn*fx;
    double a3=rz-nz*t-xi;
    D[1]=-0.5*(cs*gx*a1+sn*gx*a2+a3);
    return a1*a1+a2*a2+a3*a3;
  }
} c;

double cab(double *X, double *D){
  return c.cab(X, D);
}


struct asens{
  bool flag;
#define ANUM 11
  float mas, s[ANUM];
  vector< V<3> > ds;
  map<xppc::ikey, V<3> > cx;
  vector< float > cs;
  map<xppc::ikey, float > dx;

  asens(){
    {
      mas=1; s[0]=1;
      for(int i=1; i<ANUM; i++) s[i]=0;

      ifstream inFile("as", ifstream::in);
      if(!inFile.fail()){
	int n=0;
	float aux;
	if(inFile >> aux) mas=aux, n++;
	for(int i=0; i<ANUM; i++) if(inFile >> aux) s[i]=aux, n++;
	if(n>0) cerr<<"Loaded "<<n<<" angsens coefficients"<<endl;
	else{ cerr<<"File as did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else cerr<<"Could not open file as"<<endl;
    }

    {
      ifstream inFile("zs", ifstream::in);
      if(!inFile.fail()){
	int n;
	V<3> dir;
	while(inFile >> n >> dir[0] >> dir[1] >> dir[2]) ds.push_back(dir); 
	if(ds.size()>0) cerr<<"Loaded "<<ds.size()<<" directions"<<endl;
	else{ cerr<<"File zs did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else cerr<<"Could not open file zs"<<endl;
    }

    {
      ifstream inFile("cx", ifstream::in);
      if(!inFile.fail()){
	xppc::ikey om;
	V<3> dir;
	float r;
	while(inFile >> om.str >> om.dom >> dir[0] >> dir[1] >> dir[2] >> r) cx[om]=dir;
	if(cx.size()>0) cerr<<"Loaded "<<cx.size()<<" DOM orientations"<<endl;
	else{ cerr<<"File cx did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else cerr<<"Could not open file cx"<<endl;
    }

    {
      ifstream inFile("cs", ifstream::in);
      if(!inFile.fail()){
	int n;
	float dir;
	while(inFile >> n >> dir) cs.push_back(dir); 
	if(cs.size()>0) cerr<<"Loaded "<<cs.size()<<" directions"<<endl;
	else{ cerr<<"File cs did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else cerr<<"Could not open file cs"<<endl;
    }

    {
      ifstream inFile("dx", ifstream::in);
      if(!inFile.fail()){
	xppc::ikey om;
        float dir;
	float r;
	while(inFile >> om.str >> om.dom >> dir >> r) dx[om]=dir;
	if(dx.size()>0) cerr<<"Loaded "<<dx.size()<<" cable positions"<<endl;
	else{ cerr<<"File dx did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else cerr<<"Could not open file dx"<<endl;
    }

    flag=!ds.empty() || !cs.empty();

    if(!cx.empty() || !ds.empty() || !dx.empty() || !cs.empty()) cerr<<"Warning: local DOM effects might be applied twice!"<<endl;
  }

  float f(float x, float z=NAN){
    float ret=1, sum=s[0];
    if(mas>0){
      {
	float y=1;
	for(int i=1; i<ANUM; i++){ y*=x; sum+=s[i]*y; }
      }

      ret=min(1.f, sum/mas);
    }
    else if(!isnan(z)){
      ret=-z>sum*OMR?1:0;
    }
    return ret;
  }

  float f(const xppc::pout & p, const xppc::ikey& om, int n = -1, int m = -1){
    const float * h = p.n, * r = p.r;
    {
      float ph=NAN;
      if(m>=0) ph=cs[m];
      else{
	map<xppc::ikey, float >::iterator di=dx.find(om);
	if(di!=dx.end()) ph = di->second;
      }

      if(!isnan(ph)){
	if(c.cab(h, r, ph)) return 0;
      }
    }

    if(n>=0){
      V<3> & d = ds[n];
      return f(h[0]*d[0]+h[1]*d[1]+h[2]*d[2], r[0]*d[0]+r[1]*d[1]+r[2]*d[2]);
    }

    {
      map<xppc::ikey, V<3> >::iterator ci=cx.find(om);
      if(ci!=cx.end()){
	V<3> & d = ci->second;
	return f(h[0]*d[0]+h[1]*d[1]+h[2]*d[2], r[0]*d[0]+r[1]*d[1]+r[2]*d[2]);
      }
      else return f(h[2], r[2]);
    }
  }
} ang;


struct hit{
  float t, q;

  hit() : t(0), q(0){ }

  hit(float t, float q){
    this->t=t; this->q=q;
  }
};

struct hix : hit{
  float c;

  hix() : hit(), c(NAN){ }

  hix(float t, float q, float c) : hit(t, q){
    this->c=c;
  }
};

typedef pair<int, int> key;
struct xkey:key{
  int pmt;

  xkey(): key(), pmt(0){}
  xkey(key x): key(x), pmt(0){}

  xkey(key x, int y): key(x){
    pmt=sumw?0:y;
  }
};

ostream& operator<<(ostream& s, const xkey& x) {
  s<<x.first<<" "<<x.second;
  if(xppc::nextgen) s<<"_"<<x.pmt;
  return s;
}

key fla;
float fln[3]={0,0,1};

typedef struct{
  double r[3];
  xkey dom;

  double t0;
  vector< pair<double, double> > ert;

  bool isok(double t){
    for(vector< pair<double, double> >::const_iterator i=ert.begin(); i!=ert.end(); ++i) if(i->first<=t && t<i->second) return false;
    return true;
  }

  vector<hit> draw, sraw;
  double dbin[num], sbin[num], dtot, stot, qtot;

  int n1, n2, ntop;
  double *dtop, *stop;

  void prep(double t0){
    this->t0=t0;
    for(int k=0; k<num; k++) dbin[k]=0; dtot=0; qtot=0;
    for(vector<hit>::const_iterator j=draw.begin(); j!=draw.end(); ++j){
      qtot+=j->q;
      if(isok(j->t)){
	double t=floor((j->t-t0)/lbin), q=j->q;
	if(cutz) if(t<=0 || t>=num-2) q=0;
	if(t<0) t=0; else if(t>=num) t=num-1;
	dbin[int(t)]+=q; dtot+=q;
      }
    }

    {
      int dmax=0;
      double vmax=0;

      for(int i=0; i<num; i++) if(dbin[i]>vmax) vmax=dbin[i], dmax=i;
      n1=max(0, dmax-(int) round(500/lbin));

      for(int i=num-1; i>dmax; i--) if(dbin[i]==vmax){ dmax=i; break; }
      n2=min(dmax+(int) round(1000/lbin), num);

      double d1=0, d2=0;
      for(int i=0; i<n1; i++) d1+=dbin[i];
      for(int i=n2; i<num; i++) d2+=dbin[i];
      dbin[n1]+=d1; dbin[n2-1]+=d2;

      dtop=&dbin[n1], stop=&sbin[n1], ntop=n2-n1;
    }

    optimize();
  }

  void fill(double dt){
    if(isfinite(dt)){
      for(int k=n1; k<n2; k++) sbin[k]=0; stot=0;
      for(vector<hit>::const_iterator j=sraw.begin(); j!=sraw.end(); ++j) if(isok(j->t+dt+t0)){
	double t=floor((j->t+dt)/lbin), q=j->q;
	if(cutz) if(t<=0 || t>=num-2) q=0;
	if(t<n1) t=n1; else if(t>=n2) t=n2-1;
	sbin[int(t)]+=q; stot+=q;
      }

      join();
    }
  }

  double llh(double dt, double nsim){
    fill(dt);
    double result=0;
    if(norm){
      if(stot>0) nsim=stot*drep/dtot;
      else return 0;
    }
    for(int k=0; k<fnum; k++) result+=llf(sim[k], dat[k], nsim, drep);
    return result;
  }

  void clear(){
    sraw.clear();
  }

  double G(double L, double D){
    const int prior=0;  // 0 or 8
    double R=D<=0?0:D*log(D/L);
    double Q=prior>0?-log(prior):gsl_sf_lnfact(L);
    return R+Q;
  }

  double opt[num];
  int lastchange[num];

  double dat[num], sim[num];
  int fnum;

  void optimize(){
    for(int n=0; n<ntop; n++){
      double sD=0, mxo;
      int mxj=-1;
      for(int j=n; j>=0; j--){
	sD+=dtop[j];
	int len=n-j+1;
	double aux=G(len, sD);
	if(j>0) aux+=opt[j-1];
	if(mxj==-1 || aux>mxo){ mxo=aux; mxj=j; }
      }
      opt[n]=mxo; lastchange[n]=mxj;
    }

    int f=ntop-1; fnum=0;
    for(int k=lastchange[f]; ; k=lastchange[k-1]){
      double d=0;
      for(; f>=k; f--) d+=dtop[f];
      dat[fnum++]=d; if(k==0) break;
    }
  }

  void join(){
    int f=ntop-1, fnum=0;
    for(int k=lastchange[f]; ; k=lastchange[k-1]){
      double s=0;
      for(; f>=k; f--) s+=stop[f];
      sim[fnum++]=s; if(k==0) break;
    }
  }

  void bins(){
    int f=ntop-1;
    double q=0;
    for(int k=lastchange[f]; ; k=lastchange[k-1]){
      double s=0;
      if(f==ntop-1) s+=num-n2;
      for(; f>=k; f--) s++;
      if(f==0) s+=n1;
      noise.push_back(s), q+=s;
      if(k==0) break;
    }
    noise_q.push_back(q);
  }

  static double llf(double s, double d, double nsim, double ndat, double sig = LSSL::SIG){
    if(s<0 || d<0){ cerr<<"Error!"<<endl; exit(1); }
    if(s==0 && d==0) return 0;

    double s2=sig*sig, tol=1.e-8;
    double sd=s+d, ntot=ndat+nsim;
    double ls=sd/ntot, ld=ls, llx;

    if(sig>0){
      double dls;
      do{
	double f=nsim*ls-s-log(ld/ls)/s2;
	double df=1+(1/(ndat*ld)+1/(nsim*ls))/s2;
	dls=-f/(nsim*df); ls+=dls; ld=(sd-nsim*ls)/ndat;
	if(ls<0 || ld<0){ cerr<<"Error!"<<endl; exit(1); }
      } while(fabs(dls)>tol*ls);

      llx=square(log(ld/ls))/(2*s2);
    }
    else llx=0;

    double lls=0; if(s>0) lls+=s*log((s/nsim)/ls);
    double lld=0; if(d>0) lld+=d*log((d/ndat)/ld);

    return llx+lls+lld;
  }

  bool in(){
    return qtot<qsat*drep && qtot>qmin*drep && (!flasher || dom.first!=abs(fla.first) || abs(dom.second-fla.second)>fsep);
  }
} dat;

set<key> bad;
bool isok(const key & om){
  return bad.find(om)==bad.end();
}

map<xkey, dat> all;
unsigned int nchn=1, ndof=1;

double COG[3]={0};
double qdat=0, qsim=0, nohit=0, tshift=0;
bool unh=false;

double qint(double q){
  // q=round(q); if(q<=0) q=1;
  return q;
}

int lcsp=0; // 2
bool spef=false, nois=false;

float spe(){
  if(spef){
    float q=-1.f, pe=0.2987f, s0=0.2916f, qt=0.5057f;
    while(!isfinite(q) || q<0.f || q>3.f) q=xppc::xrnd()<pe?-qt*logf(1.f-xppc::xrnd()):1.f+s0*xppc::grnd();
    return q;
  }
  else return 1.;
}

float spe(unsigned int n){
  if(spef){
    double q=0;
    for(unsigned int i=0; i<n; i++) q+=spe();
    return q;
  }
  else return n;
}

float spe(const pair<xppc::ihit, vector<xppc::pout> >& dh, int n = -1, int m = -1){
  unsigned int q=0;
  const vector<xppc::pout>& hits = dh.second;
  const xppc::ikey& om = dh.first.omkey;
  for(vector<xppc::pout>::const_iterator i=hits.begin(); i!=hits.end(); ++i) if(xrnd()<ang.f(*i, om, n, m)) q++;
  return spe(q);
}

float delay(){
  const double width=fdur;
  double delay=width*xrnd();
  if(oldf) if(xrnd()<exp(-delay/tau)) delay+=width; // for tau << fdur only!
  return delay;
}

struct xyz{
  double x, y, z;

  xyz(){ x=0, y=0, z=0; }

  xyz(double x, double y, double z){
    this->x=x, this->y=y, this->z=z;
  }

  void ini(){
    {
      FILE * f = fopen("bad", "r");
      if(f){
	int str, dom;
	while(fscanf(f, "%d %d", &str, &dom)==2) bad.insert(make_pair(str, dom));
	fclose(f);
	cerr<<"Bad DOM list read; no-hit DOMs treatment enabled"<<endl;
	unh=true;
      }
    }

    {
      FILE * f = fopen("ert", "r");
      if(f){
	int str, dom;
	double ti, tf;
	while(fscanf(f, "%d %d %lf %lf", &str, &dom, &ti, &tf)==4) all[make_pair(str, dom)].ert.push_back(make_pair(ti, tf));
	fclose(f);
	cerr<<"Errata time interval list read; masking bad pulse data"<<endl;
      }
    }

    {
      const double dt=0;
      double t0=NAN, tx=NAN;

      FILE * f = fopen("dat", "r");
      if(f){
	int str, dom, pmt=0;
	double t, q;
	while(xppc::nextgen ?
	      fscanf(f, "%d %d_%d %lf %lf", &str, &dom, &pmt, &t, &q)==5 :
	      fscanf(f, "%d %d %lf %lf", &str, &dom, &t, &q)==4){
	  xkey om(make_pair(str, dom), pmt);
	  xppc::ikey ik; ik.str=str; ik.dom=dom;
	  if(ik.isinice() && isok(om)){
	    if(std::isnan(t0) || t<t0) t0=t;
	    if(std::isnan(tx) || t>tx) tx=t;
	    q=qint(q*drep);
	    if(q>0) all[om].draw.push_back(hit(t, q));
	  }
	}
	fclose(f);
      }

      {
	double aux=(tx-t0+1)/num; lbin=aux>bin?aux:bin;
	cerr<<"T="<<(tx-t0)<<" "<<(num*bin)<<" lbin="<<lbin<<" ";
      }

      tshift=t0-dt; nchn=1, ndof=1;
      for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	dat & d = i->second; d.dom=i->first;
	d.prep(tshift); if(d.in()) nchn++, ndof+=d.fnum;
      }

      unsigned int nhit=0;
      for(vector<xppc::OM>::const_iterator i=xppc::i3oms.begin(); i!=xppc::i3oms.end(); ++i){
	key om=make_pair(i->str, i->dom);
	if(i->isinice() && isok(om) && all.find(om)==all.end()){
	  float h;
	  if(xppc::hvs.empty()) h=1200;
	  else{
	    map<xppc::ikey, float>::iterator j=xppc::hvs.find(*i);
	    h=j==xppc::hvs.end()?0:j->second;
	  }

	  if(h>0) nhit++;
	}
      }

      noise.push_back(nhit*num);
      noise_q.push_back(noise.front());
      for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	dat & d=i->second;
	if(d.in()) d.bins();
      }

      nsq=nsrt*lbin*drep;
      for(vector<double>::iterator i=noise.begin(); i!=noise.end(); ++i) *i*=nsq;
      for(vector<double>::iterator i=noise_q.begin(); i!=noise_q.end(); ++i) *i*=nsq;
    }

    {
      map<key, xyz> geo;

      for(vector<xppc::OM>::const_iterator i=xppc::i3oms.begin(); i!=xppc::i3oms.end(); ++i)
	geo.insert(make_pair(make_pair(i->str, i->dom), xyz(i->r[0], i->r[1], i->r[2])));

      x=0, y=0, z=0; qdat=0;
      for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	xyz g = geo[i->first];
	dat & d=i->second;
	double q=d.qtot;
	d.r[0]=g.x, d.r[1]=g.y, d.r[2]=g.z;
	x+=q*g.x; y+=q*g.y; z+=q*g.z; qdat+=q;
      }
      if(qdat>0) x/=qdat, y/=qdat, z/=qdat;
    }

    COG[0]=x, COG[1]=y, COG[2]=z;
    cerr<<"COG="<<x<<","<<y<<","<<z<<endl<<"Ndof="<<ndof<<" Qdat="<<qdat<<endl;
  }
};

struct dir{
  double th, ph, n[3];
  double p1[3], p2[3], u[3];

  dir() : th(0), ph(0){ n[0]=0; n[1]=0; n[2]=1; }

  void normalize(){
    double r=0;  for(int i=0; i<3; i++) r+=n[i]*n[i];
    r=1/sqrt(r); for(int i=0; i<3; i++) n[i]*=r;
  }

  double xvmf(double da){
    double ka=square(da*cv), xi;
    do{ xi=1+ka*log(xrnd()); } while (xi<-1);
    return acos(xi);
  }

  void setn(){
    double th=cv*this->th, ph=cv*this->ph;
    double st=sin(th), ct=cos(th);
    double sp=sin(ph), cp=cos(ph);
    n[0]=-st*cp, n[1]=-st*sp, n[2]=-ct;
  }

  void seta(){
    th=acos(-n[2])/cv; ph=atan2(-n[1], -n[0])/cv;
  }

  void setp(){
    {
      int i0=0;
      double n0=n[0];
      for(int i=1; i<3; i++) if(fabs(n[i])>fabs(n0)){ i0=i; n0=n[i]; }
      int i1=(i0+1)%3, i2=(i0+2)%3;
      double n1=n[i1], n2=n[i2], nq=n0*n0;
      double r1=1/sqrt(nq+n1*n1), r2=1/sqrt(nq+n2*n2);

      p1[i0]=-n1*r1; p1[i1]=n0*r1; p1[i2]=0;
      p2[i0]=-n2*r2; p2[i1]=0; p2[i2]=n0*r2;
    }
    {
      double q1[3], q2[3], r1=0, r2=0;
      for(int i=0; i<3; i++){
	double a1=p1[i]-p2[i]; q1[i]=a1; r1+=a1*a1;
	double a2=p1[i]+p2[i]; q2[i]=a2; r2+=a2*a2;
      }

      r1=1/sqrt(r1); r2=1/sqrt(r2);

      for(int i=0; i<3; i++){
	p1[i]=q1[i]*r1;
	p2[i]=q2[i]*r2;
      }
    }
  }

  void setu(){
    double zi=xrnd(); zi*=2*M_PI;
    double px=cos(zi), py=sin(zi);
    for(int i=0; i<3; i++) u[i]=px*p1[i]+py*p2[i];
  }

  void rotate(double z){
    if(flasher) return;
    double cs=cos(z), si=sin(z);
    {
      double r=0;
      for(int i=0; i<3; i++){
	double a=n[i]*cs+u[i]*si;
	n[i]=a; r+=a*a;
      }
      {
	r=1/sqrt(r);
	for(int i=0; i<3; i++) n[i]*=r;
      }
    }
  }
};

struct cascade : xyz, dir{
  double e, t, s, a;
  vector<double> unfE;

  cascade() : xyz(), dir(){
    e=0, t=0, s=1, a=1;
  }

  bool ini(){
    xyz::ini();
    bool flag=false;

    ifstream h("ini", ifstream::in);
    if(!h.fail()){
      string in;
      if(getline(h, in)){
	cascade c;
	int num=sscanf(in.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &c.x, &c.y, &c.z, &c.th, &c.ph, &c.e, &c.t, &c.s, &c.a);
	if(num>=5){
	  c.t-=tshift; c.setn(); * this = c;
	  if(!flasher || x!=0 || y!=0 || z!=0) flag=true;
	  if(mlpd && !flasher){
	    double dx=COG[0]-x, dy=COG[1]-y, dz=COG[2]-z;
	    double dot=dx*n[0]+dy*n[1]+dz*n[2];
	    x+=dot*n[0], y+=dot*n[1], z+=dot*n[2];
	    if(num>=7) t+=dot/sol; else t=-tshift;
	  }
	  cerr<<"Setting initial cascade parameters from file ini."<<endl;
	}
      }
      if(getline(h, in)){
	stringstream str(in);
	double aux;
	while(str>>aux) unfE.push_back(aux);
	if(!unfE.empty()) cerr<<"Setting unfolded pattern with "<<unfE.size()<<" entries from file ini."<<endl;
      }
      if(getline(h, in)){
	double r, dr, da, di, lx;
	int num=sscanf(in.c_str(), "%lf %lf %lf %lf %lf", &r, &dr, &da, &di, &lx);
	if(num>=4){
	  rr=r, dR=dr, dA=da, dI=di;
	  cerr<<"Setting search parameters from file ini."<<endl;
	}
	if(num>=5){
	  LX=lx;
	  cerr<<"Setting LX="<<LX<<endl;
	}
      }
      h.close();
    }

    else{
      char * ang=getenv("ANGR");
      if(ang!=NULL){
	string in(ang);
	double nx, ny, nz;
	int num=sscanf(in.c_str(), "%lf %lf %lf", &nx, &ny, &nz);
	if(num==3){
	  n[0]=nx, n[1]=ny, n[2]=nz, dA=0; seta();
	  cerr<<"Setting ANGR="<<in<<endl;
	}
      }
    }

    step();
    return flag;
  }

  double gaus(){
    static bool flag=false;
    static double r1, r2;
    flag=!flag;
    if(flag){
      double ra=sqrt(-2*log(xrnd())), rb=2*M_PI*xrnd();
      r1=ra*cos(rb); r2=ra*sin(rb);
      return r1;
    }
    else{
      flag=false;
      return r2;
    }
  }

  void move(double dr=dR, double da=dA, double di=dI){
    x+=rr*n[0], y+=rr*n[1], z+=rr*n[2];

    dr/=M_SQRT3;
    x+=dr*gaus(), y+=dr*gaus(), z+=dr*gaus();

    setu(); rotate(xvmf(da));

    x-=rr*n[0], y-=rr*n[1], z-=rr*n[2];

    s*=exp(di*gaus()), a*=exp(di*gaus());
  }

  void unfold(double[], double[], vector< pair<xkey, hit> >[], int);
  void unfold_t(double[], double[], vector< pair<xkey, hit> >[], int, double);

  void runc(map< xkey, vector<hix> > &);
  void runf(map< xkey, vector<hix> > &);
  void ffit();
  double optte(double&, double&);
  void ppc(vector< pair<xkey, hix> > *);
  void end();

  double llh(vector< pair<xkey, hix> > *, bool, bool);
  void printout(vector< pair<xkey, hix> > [], int);
  void refine(const multimap<double, cascade> &, int);
};

ostream& operator<<(ostream& o, cascade& c){
  c.seta();
  return o<<c.x<<" "<<c.y<<" "<<c.z<<" "<<c.th<<" "<<c.ph<<" "<<c.e<<" "<<c.t+tshift<<" "<<c.s<<" "<<c.a;
}

double fllh(double dt, double nsim){
  double llh=norm?0:dat::llf(nohit, 0, nsim, drep);
  for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
    dat & d=i->second;
    if(d.in()) llh+=d.llh(dt, nsim);
  }
  return llh/ndof;
}

double fllh(const gsl_vector *x, void *p){
  return fllh(gsl_vector_get(x, 0), srep*exp(gsl_vector_get(x, 1)));
}

double fllt(double x, void * p)
{
  (void)(p);
  return fllh(x, srep);
}

void cascade::unfold(double Q[], double X[], vector< pair<xkey, hit> > sto[], int jnum){
  double * A = new double[nchn*jnum]();

  for(int m=0; m<jnum; m++){
    vector< pair<xkey, hit> > & s = sto[m];
    for(vector< pair<xkey, hit> >::iterator h=s.begin(); h!=s.end(); ++h){
      const hit& tq = h->second;
      double t=tq.t, q=tq.q;
      map<xkey, dat>::iterator i=all.find(h->first);
      if(i==all.end()){ if(unh) nohit+=q; }
      else i->second.sraw.push_back(hit(t, q));
    }
    {
      int n=0;
      A[n+++m*nchn]=nohit;
      for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	dat & d=i->second;
	if(d.in()){
	  d.fill(0);
	  double q=0;
	  for(int k=0; k<d.fnum; k++) q+=d.sim[k];
	  A[n+++m*nchn]=q;
	}
      }
      end();
    }
  }

  {
    // NNLS::nnls(A, Q, X, nchn, jnum);
    LSSL::wref(A, Q, X, noise_q.data(), nchn, jnum, -2);
  }

  delete A;
}

void cascade::unfold_t(double D[], double X[], vector< pair<xkey, hit> > sto[], int jnum, double dt=0){
  bool time=false, same=false, done=false;
  double * A = new double[ndof*jnum]();
  double del=4*lbin, low=lbin/25, best=NAN;
  double sdt=dt;

  while(true){
    for(int m=0; m<jnum; m++){
      vector< pair<xkey, hit> > & s = sto[m];
      for(vector< pair<xkey, hit> >::iterator h=s.begin(); h!=s.end(); ++h){
	const hit& tq = h->second;
	double t=tq.t, q=tq.q;
	map<xkey, dat>::iterator i=all.find(h->first);
	if(i==all.end()){ if(unh) nohit+=q; }
	else i->second.sraw.push_back(hit(t, q));
      }
      {
	int n=0;
	A[n+++m*ndof]=nohit;
	for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	  dat & d=i->second;
	  if(d.in()){
	    d.fill(dt);
	    for(int k=0; k<d.fnum; k++) A[n+++m*ndof]=d.sim[k];
	  }
	}
	end();
      }
    }

    if(time){
      double rnorm=NNLS::nnls(A, D, X, ndof, jnum);
      // double rnorm=NNLS::nsum(A, D, X, ndof, jnum);
      // double rnorm=LSSL::wref(A, D, X, noise.data(), ndof, jnum, -2);
      // double rnorm=LSSL::wllh(A, D, X, noise.data(), ndof, jnum);
      // cout<<dt<<" "<<rnorm<<" "<<del<<endl;

      if(done) break;

      if(std::isnan(best)) best=rnorm;
      else{
	if(fabs(del)<low){
	  if(rnorm<best){
	    best=rnorm;
	    break;
	  }
	  else{
	    del=sdt-dt;
	    done=true;
	  }
	}
	else{
	  if(rnorm<best){
	    best=rnorm; sdt=dt;
	    if(same) del*=2; else same=true;
	  }
	  else{
	    dt-=del; del/=-2;
	    if(same){ del/=2; same=false; }
	  }
	}
      }
      dt+=del;
    }
    else break;
  }

  {
    // NNLS::nnls(A, D, X, ndof, jnum);
    LSSL::wref(A, D, X, noise.data(), ndof, jnum, -2);
  }

  delete A;
  if(time) t+=dt;
}

void cascade::runc(map< xkey, vector<hix> > & sim){
  pair<int, unsigned long long> id=make_pair(0, 0);
  const double ini=1.e6, emax=1.e8;

  if(oldf){
    if(e<=0) e=ini;
    xppc::sett(n[0], n[1], n[2], id, 0);
    xppc::addp(x, y, z, t, e, 0, 10);
    xppc::eout();

    for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
      const xppc::ihit & h = i->first;
      xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
      sim[om].push_back(hix(h.time, spe(*i), h.dir));
    }

    xppc::efin();
  }
  else{
    if(mlpd){
      double dx=COG[0]-x, dy=COG[1]-y, dz=COG[2]-z;
      double dot=dx*n[0]+dy*n[1]+dz*n[2];
      x+=dot*n[0], y+=dot*n[1], z+=dot*n[2], t+=dot/sol;
    }

    const double etot=ini*srep; // total simulated energy
    const double step=7.5;  // distance between simulated cascades
    const double dtrs=100;  // max. distance to the closest DOM at ends
    const double dpad=100;  // padding to the segment of track with hits

    double dmin=0, dmax=0;  // fit end points

    if(mlpd) for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
      dat & d=i->second;

      double dx=d.r[0]-x, dy=d.r[1]-y, dz=d.r[2]-z;
      double rn=dx*n[0]+dy*n[1]+dz*n[2];
      dx-=rn*n[0], dy-=rn*n[1], dz-=rn*n[2];
      double dd=sqrt(dx*dx+dy*dy+dz*dz);

      if(dd<dtrs){
	if(rn<dmin) dmin=rn;
	else if(rn>dmax) dmax=rn;
      }
    }

    int jmin=(int) round((dmin-dpad)/step);
    int jmax=(int) round((dmax+dpad)/step);
    int jnum=jmax-jmin+1;  // cout<<jnum<<" "<<ndof<<endl;

    if(!mlpd) jmin=0, jmax=0, jnum=1;

    if(unfE.empty()){
      double D[ndof], Q[nchn], E[jnum], X[jnum], T[jnum];

      {
	int n=0, m=0;
	D[n++]=0; Q[m++]=0;
	for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	  dat & d=i->second;
	  if(d.in()){
	    double q=0;
	    for(int k=0; k<d.fnum; k++) D[n++]=d.dat[k], q+=d.dat[k];
	    Q[m++]=q;
	  }
	}

	double e0=etot/jnum;
	e=0; for(int j=0; j<jnum; j++) E[j]=e0, T[j]=E[j];
      }

      double dt=0;
      vector< pair<xkey, hit> > sto[jnum];

      for(int l=0; l<2; l++){
	// if(l==1) for(int j=0; j<jnum; j++) sto[j].clear(), T[j]=0;
	size_t qx[jnum]; for(int i=0; i<jnum; i++) qx[i]=i;
	gsl_ran_shuffle(xr.r, qx, jnum, sizeof(size_t));

	for(int J=0; J<jnum; J++){
	  int j=qx[J];

	  double rn=step*(j+jmin);
	  xppc::sett(n[0], n[1], n[2], id, j);
	  if(e>0) E[j]*=etot/e, T[j]+=E[j];
	  xppc::addp(x+rn*n[0], y+rn*n[1], z+rn*n[2], t+rn/sol, E[j], mlpd?step:0, 10);
	}
	xppc::eout();

	qsim=0;
	for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
	  const xppc::ihit & h = i->first;
	  xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
	  if(isok(om)){
	    hit tq(h.time, spe(*i));
	    sto[h.track.frame].push_back(make_pair(om, tq));
	    if(l>0 && !fast){
	      map<xkey, dat>::iterator i=all.find(om);
	      if(i==all.end()){ if(unh) nohit+=tq.q; }
	      else i->second.sraw.push_back(tq);
	      qsim+=tq.q;
	    }
	  }
	}
	xppc::efin();

	if(qsim>0){
	  double t=0, e=1;
	  double llh=optte(t, e);
	  end();
	  dt=t;
	  unfold_t(D, X, sto, jnum, t);
	}
	else unfold(Q, X, sto, jnum);

	e=0; for(int j=0; j<jnum; j++) E[j]=T[j]*X[j]/drep, e+=E[j];
	// cout<<"E:"; for(int j=0; j<jnum; j++) cout<<" "<<E[j]; cout<<endl;
	if(e>emax){
	  cerr<<"Warning: e="<<e<<" exceeded max energy of "<<emax<<"!"<<endl;
	  if(fail) exit(2);
	  for(int j=0; j<jnum; j++) E[j]*=emax/e; e=emax;
	}
      }
      t+=dt;
      for(int j=0; j<jnum; j++) unfE.push_back(e>0?E[j]/e:0);
      cout<<"E:"; for(vector<double>::iterator i=unfE.begin(); i!=unfE.end(); ++i) cout<<" "<<*i; cout<<endl;
    }

    {
      xppc::sett(n[0], n[1], n[2], id, 0);
      for(int j=0; j<jnum; j++){
	double rn=step*(j+jmin);
	xppc::addp(x+rn*n[0], y+rn*n[1], z+rn*n[2], t+rn/sol, e*unfE[j], mlpd?step:0, 10);
      }
      xppc::eout();

      for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
	const xppc::ihit & h = i->first;
	xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
	sim[om].push_back(hix(h.time, spe(*i), h.dir));
      }

      xppc::efin();
    }
  }
}

bool pass=false;

void cascade::ffit(){
  if(pass && ang.flag){
    for(int n=ang.ds.empty()?-1:0; n<(int) ang.ds.size(); n++){
      for(int m=ang.cs.empty()?-1:0; m<(int) ang.cs.size(); m++){
	map< xkey, vector<hix> > sim;
	for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
	  const xppc::ihit & h = i->first;
	  xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
	  sim[om].push_back(hix(t+h.time+delay(), spe(*i, n, m), h.dir));
	}

	for(map< xkey, vector<hix> >::const_iterator i=sim.begin(); i!=sim.end(); ++i){
	  const xkey& om = i->first;

	  if(isok(om)){
	    const vector<hix>& ts = i->second;
	    map<xkey, dat>::iterator i=all.find(om);
	    for(vector<hix>::const_iterator j=ts.begin(); j!=ts.end(); ++j){
	      const hix & h = *j;

	      if(i==all.end()){ if(unh) nohit+=h.q; }
	      else i->second.sraw.push_back(h);
	      qsim+=h.q;
	    }
	  }
	}

	for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	  const xkey & om=i->first;
	  dat & d=i->second;
	  if(d.in()) cout<<n<<" "<<m<<" "<<om<<" "<<d.llh(0, srep)/ndof<<endl;
	}
	end();
      }
    }
  }
}

void cascade::runf(map< xkey, vector<hix> > & sim){
  const double ini=5, emax=100, bunch=2.5e9;

  if(oldf){
    if(e<=0) e=ini;
    xppc::flone(e*bunch*(ang.flag?srep:1));
    xppc::eout();
    ffit();

    for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
      const xppc::ihit & h = i->first;
      xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
      sim[om].push_back(hix(t+h.time+delay(), spe(*i), h.dir));
    }

    xppc::efin();
  }
  else{
    pair<int, unsigned long long> id=make_pair(0, 0);

    const double etot = (e>0 && munf>=10 ? e : ini) * srep;  // total number of photon bunches
    const double step = 5.;  // angular distance between LED beam directions [degrees]
    int j0, j1, j2, jnum;  // angular bins

    /*  munf:
	0:1 component (mlpd=false)
	1:2+72 down/up+72 5-deg bins
	2:2+1 (down/up+led)
	3:2+1 (back/forward+led)
	4:1+1 (isotropic+led)
	+10:take e from file "ini", otherwise ini=5
    */

    switch(munf%10){
    case 1:  j0=2, j1=-1, j2=-2; jnum=j0+360/step; // (2 linear components + 5-deg bins)
      break;
    case 2:
    case 3:  j0=2, j1=-1, j2=-4; jnum=j0+1; // (2 linear components + LED beam)
      break;
    case 4:  j0=1, j1=-3, j2=-4; jnum=j0+1; // (1 isotrop component + LED beam)
      break;
    case 0:
    default: j0=0, j1=-4, j2=-4; jnum=j0+1; // one flset-defined beam
    }

    if(unfE.empty()){
      double D[ndof], Q[nchn], E[jnum], X[jnum], T[jnum];

      {
	int n=0, m=0;
	D[n++]=0; Q[m++]=0;
	for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	  dat & d=i->second;
	  if(d.in()){
	    double q=0;
	    for(int k=0; k<d.fnum; k++) D[n++]=d.dat[k], q+=d.dat[k];
	    Q[m++]=q;
	  }
	}

	double e0=etot/jnum;
	e=0; for(int j=0; j<jnum; j++) E[j]=e0, T[j]=E[j];
      }

      double dt=0;
      vector< pair<xkey, hit> > sto[jnum];

      for(int l=0; l<2; l++){
	if(mlpd){
	  size_t qx[jnum]; for(int i=0; i<jnum; i++) qx[i]=i;
	  gsl_ran_shuffle(xr.r, qx, jnum, sizeof(size_t));

	  for(int J=0; J<jnum; J++){
	    int j=qx[J];

	    float n[3], r[3]={(float) x, (float) y, (float) z};
	    if(j1==-1 && j<j0){
	      // n[0]=0, n[1]=0, n[2]=2*j-1;
	      for(int i=0; i<3; i++) n[i]=(2*j-1)*fln[i];
	    }
	    else if(j2==-2){
	      sincosf(cv*step*(j-j0), &n[1], &n[0]); n[2]=0;
	      xppc::flshift(r, n, fln);
	    }
	    else{
	      n[0]=0, n[1]=0, n[2]=0;
	    }

	    xppc::sett(n[0], n[1], n[2], id, j);
	    if(e>0) E[j]*=etot/e, T[j]+=E[j];
	    xppc::addp(r[0], r[1], r[2], 0.f, E[j]*bunch, 0.f, j<j0?j1:j2);
	  }
	}
	else{
	  for(int j=0; j<jnum; j++){
	    if(e>0) E[j]*=etot/e, T[j]+=E[j];
	    xppc::flone(E[j]*bunch);
	  }
	}
	xppc::eout();

	qsim=0;
	for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
	  const xppc::ihit & h = i->first;
	  xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
	  if(isok(om)){
	    hit tq(t+h.time+delay(), spe(*i));
	    sto[h.track.frame].push_back(make_pair(om, tq));
	    if(l>0 && !fast){
	      map<xkey, dat>::iterator i=all.find(om);
	      if(i==all.end()){ if(unh) nohit+=tq.q; }
	      else i->second.sraw.push_back(tq);
	      qsim+=tq.q;
	    }
	  }
	}
	xppc::efin();

	if(qsim>0){
	  double t=0, e=1;
	  double llh=optte(t, e);
	  end();
	  dt=t;
	  unfold_t(D, X, sto, jnum, t);
	}
	else unfold(Q, X, sto, jnum);

	e=0; for(int j=0; j<jnum; j++) E[j]=T[j]*X[j]/drep, e+=E[j];
	// cout<<"E:"; for(int j=0; j<jnum; j++) cout<<" "<<E[j]; cout<<endl;
	if(e>emax){
	  cerr<<"Warning: e="<<e<<" exceeded max energy of "<<emax<<"!"<<endl;
	  if(fail) exit(2);
	  for(int j=0; j<jnum; j++) E[j]*=emax/e; e=emax;
	}
      }
      t+=dt;
      for(int j=0; j<jnum; j++) unfE.push_back(e>0?E[j]/e:0);
      cout<<"E:"; for(vector<double>::iterator i=unfE.begin(); i!=unfE.end(); ++i) cout<<" "<<*i; cout<<endl;
    }

    {
      if(mlpd){
	for(int j=0; j<jnum; j++){
	  float n[3], r[3]={(float) x, (float) y, (float) z};
	  if(j1==-1 && j<j0){
	    // n[0]=0, n[1]=0, n[2]=2*j-1;
	    for(int i=0; i<3; i++) n[i]=(2*j-1)*fln[i];
	  }
	  else if(j2==-2){
	    sincosf(cv*step*(j-j0), &n[1], &n[0]); n[2]=0;
	    xppc::flshift(r, n, fln);
	  }
	  else{
	    n[0]=0, n[1]=0, n[2]=0;
	  }

	  xppc::sett(n[0], n[1], n[2], id, j);
	  xppc::addp(r[0], r[1], r[2], 0.f, e*unfE[j]*bunch*(ang.flag?srep:1), 0.f, j<j0?j1:j2);
	}
      }
      else{
	for(int j=0; j<jnum; j++){
	  xppc::flone(e*unfE[j]*bunch*(ang.flag?srep:1));
	}
      }
      xppc::eout();

      ffit();

      for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
	const xppc::ihit & h = i->first;
	xkey om(make_pair(h.omkey.str, h.omkey.dom), h.pmt);
	sim[om].push_back(hix(t+h.time+delay(), spe(*i), h.dir));
      }

      xppc::efin();
    }
  }
}

void cascade::ppc(vector< pair<xkey, hix> > * sav = NULL){
  qsim=0;

  // nois, spef, and lcsp are not supported with multi-PMT sensors

  for(int i=0; i<(ang.flag?1:srep); i++){
    map< xkey, vector<hix> > sim; flasher?runf(sim):runc(sim);

    if(nois){
      for(vector<xppc::OM>::const_iterator i=xppc::i3oms.begin(); i!=xppc::i3oms.end(); ++i){
	xkey om=make_pair(i->str, i->dom);
	if(isok(om)){
	  float h;
	  if(xppc::hvs.empty()) h=1200;
	  else{
	    map<xppc::ikey, float>::iterator j=xppc::hvs.find(*i);
	    h=j==xppc::hvs.end()?0:j->second;
	  }

	  if(h>0){
	    float t=-1000;
	    while((t-=log(xrnd())/nsrt)<6000) sim[om].push_back(hix(t, spe(), NAN));
	  }
	}
      }
    }

    if(spef){
      const double dthr=0.25;
      for(map< xkey, vector<hix> >::iterator i=sim.begin(), j; i!=sim.end(); ){
	const vector<hix>& qs = i->second;
	double qi=0;
	for(vector<hix>::const_iterator k=qs.begin(); k!=qs.end(); ++k) qi+=k->q;
	j=i++; if(qi<dthr) sim.erase(j);
      }
    }

    for(map< xkey, vector<hix> >::const_iterator i=sim.begin(); i!=sim.end(); ++i){
      const xkey& om = i->first;

      if(lcsp>0){
	bool skip=true;
	if(skip){
	  map< xkey, vector<hix> >::const_iterator j(i); j++;
	  if(j!=sim.end() && j->first.first==om.first && j->first.second-om.second<=lcsp) skip=false;
	}
	if(skip){
	  map< xkey, vector<hix> >::const_reverse_iterator j(i);
	  if(j!=sim.rend() && j->first.first==om.first && om.second-j->first.second<=lcsp) skip=false;
	}

	if(skip) continue;
      }

      if(isok(om)){
	const vector<hix>& ts = i->second;
	map<xkey, dat>::iterator i=all.find(om);
	for(vector<hix>::const_iterator j=ts.begin(); j!=ts.end(); ++j){
	  const hix & h = *j;

	  if(i==all.end()){ if(unh) nohit+=h.q; }
	  else i->second.sraw.push_back(h);

	  if(sav!=NULL)
	    if(isfinite(h.c)) sav[min(cnum-1, max(0, (int) floor((1+h.c)/cbin)))].push_back(make_pair(om, h));

	  qsim+=h.q;
	}
      }
    }
  }
}

void cascade::end(){
  nohit=0;
  for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i) i->second.clear();
}

double cascade::optte(double &t, double &e){
  double llh;

  if(qsim>0){
    if(norm){
      const gsl_min_fminimizer_type *T;
      gsl_min_fminimizer *s;
      gsl_function F;

      F.function = &fllt;
      F.params = 0;

      T = gsl_min_fminimizer_brent;
      s = gsl_min_fminimizer_alloc (T);

      double d = 0, a = -100, b = 100;
      gsl_min_fminimizer_set (s, &F, d, a, b);

      for(int i=0; i<100; i++){
	if(GSL_SUCCESS!=gsl_min_fminimizer_iterate(s)) break;

	d = gsl_min_fminimizer_x_minimum(s);
	a = gsl_min_fminimizer_x_lower(s);
	b = gsl_min_fminimizer_x_upper(s);

	if(GSL_SUCCESS==gsl_min_test_interval(a, b, 0.1, 0.0)) break;
      }

      gsl_min_fminimizer_free (s);

      t+=d; llh=fllh(d, srep);
      cout<<"LLH="<<llh<<" DT="<<d<<endl;
    }
    else{
    bool verbose=false;

    const int num=2;
    double xi[num]={0, 1}, ei[num]={100, 1};

    gsl_vector * x = gsl_vector_alloc(num);
    gsl_vector * u = gsl_vector_alloc(num);

    for(int i=0; i<num; i++) gsl_vector_set(x, i, xi[i]), gsl_vector_set(u, i, ei[i]);

    gsl_multimin_function F;
    F.n = num;
    F.f = fllh;
    F.params = NULL;

    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc(T, num);
    gsl_multimin_fminimizer_set(s, &F, x, u);

    int iter=0, max_iter=200;

    for(; iter<max_iter; iter++){
      if(GSL_SUCCESS!=gsl_multimin_fminimizer_iterate(s)){
	if(verbose) printf("Failed to iterate\n");
	break;
      }

      double size=gsl_multimin_fminimizer_size(s);
      if(verbose){
	cout<<iter<<" "<<s->fval<<" "<<size;
	for(int i=0; i<num; i++) cout<<" "<<gsl_vector_get(s->x, i); cout<<endl;
      }

      if(GSL_SUCCESS==gsl_multimin_test_size(size, 0.001)){
	if(verbose) printf("Converged!\n");
	break;
      }
    }

    double dt=gsl_vector_get(s->x, 0), nsim=exp(gsl_vector_get(s->x, 1));
    if(mlpd){
      const double nmax=100;
      if(nsim*nmax<1 || nsim>nmax){
	cerr<<"Warning: fitted nsim="<<nsim<<"!"<<endl;
	if(fail) exit(2);
	nsim=nsim<1?1/nmax:nmax;
      }
    }
    t+=dt; e/=nsim; llh=fllh(dt, srep);
    if(verbose) cout<<"fval="<<s->fval<<" "<<llh<<" "<<fllh(0, srep)<<endl;
    cout<<"LLH="<<llh<<" DT="<<dt<<" NS="<<nsim<<endl;

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(u);
    gsl_vector_free(x);
  }
  }
  else llh=fllh(0, srep);
  return llh;
}


double cascade::llh(vector< pair<xkey, hix> > * sav = NULL, bool first = true, bool second = true){
  double llh;

  if(first){
    pass=false;
    ppc();
    llh=optte(t, e);
    end();
  }

  if(second){
    pass=true;
    ppc(sav);
    llh=fllh(0, srep);
    end();
  }

  if(!oldf) unfE.clear();
  return llh;
}

void cascade::printout(vector< pair<xkey, hix> > sto[], int jnum){
  double * A = new double[ndof*jnum]();

  for(int m=0; m<jnum; m++){
    vector< pair<xkey, hix> > & s = sto[m];
    for(vector< pair<xkey, hix> >::iterator j=s.begin(); j!=s.end(); ++j){
      // const xkey & om=j->first;
      const hix& h = j->second;
      // cout<<"SREP "<<om<<" "<<h.t+tshift<<" "<<h.q<<" "<<srep<<endl;
      map<xkey, dat>::iterator i=all.find(j->first);
      if(i==all.end()){ if(unh) nohit+=h.q; }
      else i->second.sraw.push_back(h);
    }
    {
      int n=0;
      A[n+++m*ndof]=nohit;
      if(jnum==1){
	double llh=norm?0:dat::llf(nohit, 0, srep, drep);
	cout<<"N "<<ndof<<" "<<nohit<<" "<<llh<<endl;
      }
      for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
	dat & d=i->second;
	if(d.in()){
	  d.fill(0);
	  for(int k=0; k<d.fnum; k++) A[n+++m*ndof]=d.sim[k];
	  if(jnum==1){
	    const xkey & om=i->first;
	    double llh=d.llh(0, srep);
	    cout<<"L "<<om<<" "<<llh<<endl;
	  }
	}
      }
      end();
    }
  }

  int n=1;
  for(map<xkey, dat>::iterator i=all.begin(); i!=all.end(); ++i){
    dat & d=i->second;
    if(d.in()){
      const xkey & om=i->first;
      for(int k=0; k<d.fnum; k++, n++){
	cout<<om<<" "<<noise[n]/nsq<<" "<<d.dat[k]/drep;
	for(int m=0; m<jnum; m++) cout<<" "<<A[n+m*ndof]/srep;
	cout<<endl;
      }
    }
  }

  delete A;
}

void cascade::refine(const multimap<double, cascade> & v, int top){
  vector< V<7> > vs;

  int num=0;
  for(multimap<double, cascade>::const_iterator j=v.begin(); j!=v.end() && num<top; ++j, num++){
    const cascade & cj = j->second;
    V<7> w;
    w[0]=n[0]*cj.p1[0]+n[1]*cj.p1[1]+n[2]*cj.p1[2];
    w[1]=n[0]*cj.p2[0]+n[1]*cj.p2[1]+n[2]*cj.p2[2];
    w[2]=cj.x-x, w[3]=cj.y-y, w[4]=cj.z-z;
    w[5]=1-(n[0]*cj.n[0]+n[1]*cj.n[1]+n[2]*cj.n[2]);
    w[6]=j->first;
    vs.push_back(w);
  }

  const bool verbose=false;
  const int dim=5;
  const int siz=(dim+1)*(dim+2)/2;

  double q[num][siz], f[num];
  double par[dim][dim], val[dim];

  for(int m=0; m<num; m++){
    int k=0;
    for(int i=0; i<dim; i++){
      double xi=vs[m][i];
      q[m][k++]=xi*xi;
      for(int j=i+1; j<dim; j++, k++) q[m][k]=2*xi*vs[m][j];
      q[m][k++]=xi;
    }
    q[m][k++]=1;
    q[m][dim]-=q[m][0]; q[m][0]=vs[m][5];
    f[m]=vs[m][6];
  }

  {
    gsl_matrix * A = gsl_matrix_alloc(num, siz);
    gsl_vector * B = gsl_vector_alloc(num);

    for(int i=0; i<num; i++){
      for(int j=0; j<siz; j++) gsl_matrix_set(A, i, j, q[i][j]);
      gsl_vector_set(B, i, f[i]);
    }

    gsl_multifit_linear_workspace * W = gsl_multifit_linear_alloc(num, siz);
    gsl_vector * X = gsl_vector_alloc(siz);
    gsl_matrix * C = gsl_matrix_alloc(siz, siz);
    double c;

    gsl_multifit_linear(A, B, X, C, &c, W);

    double p[siz];
    for(int i=0; i<siz; i++) p[i]=gsl_vector_get(X, i);

    for(int i=0, k=0; i<dim; i++){
      par[i][i]=p[k++];
      for(int j=i+1; j<dim; j++, k++) par[i][j]=p[k], par[j][i]=p[k];
      val[i]=p[k++];
    }
    double a=par[0][0]/2, b=par[1][1];
    par[0][0]=a-b; par[1][1]=a+b;

    if(verbose){
      for(int i=0; i<dim; i++){ cout<<val[i]; for(int j=0; j<dim; j++) cout<<"\t"<<par[i][j]; cout<<endl; }
      cout<<p[siz-1]<<endl;
    }

    gsl_vector_free(X);
    gsl_matrix_free(C);
    gsl_multifit_linear_free(W);

    gsl_vector_free(B);
    gsl_matrix_free(A);
  }

  {
    gsl_matrix * A = gsl_matrix_alloc(dim, dim);
    gsl_vector * B = gsl_vector_alloc(dim);

    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++) gsl_matrix_set(A, i, j, 2*par[i][j]);
      gsl_vector_set(B, i, -val[i]);
    }

    gsl_vector * X = gsl_vector_alloc(dim);

    gsl_linalg_HH_solve(A, B, X);

    gsl_matrix_free(A);
    gsl_vector_free(B);

    double p[dim];
    for(int i=0; i<dim; i++) p[i]=gsl_vector_get(X, i);
    gsl_vector_free(X);

    if(verbose){ cout<<"X:"; for(int i=0; i<dim; i++) cout<<"\t"<<p[i]; cout<<endl; }

    {
      double nx=p[0], ny=p[1];
      double nz=sqrt(nx*nx+ny*ny), zx=sin(dA*cv);
      if(nz>zx) nx*=zx/nz, ny*=zx/nz, nz=zx;
      nz=sqrt(1-nz*nz);

      for(int i=0; i<3; i++) n[i]=nx*p1[i]+ny*p2[i]+nz*n[i];
    }

    {
      double dx=p[2], dy=p[3], dz=p[4];
      double dr=sqrt(dx*dx+dy*dy+dz*dz);
      if(dr>dR) dx*=dR/dr, dy*=dR/dr, dz*=dR/dr;

      x+=dx, y+=dy, z+=dz;
    }
  }
}

void loc_rs(cascade & ci){  // localized random search
  const int tries=25, steps=loop/tries, upto=steps/2;

  multimap<double, cascade> bag;
  for(int k=0; k<steps; k++){
    ci.setp();
    cascade cf;
    double best;

    multimap<double, cascade> v;
    for(int i=0; i<tries; i++){
      cascade ct=ci; ct.move();
      double temp=ct.llh();
      cout<<"|\t"<<temp<<" "<<ct<<endl;
      v.insert(make_pair(temp, ct));
      if(i==0 || temp<best) best=temp, cf=ct;
    }
    if(!mlpd && !flasher){
      cascade ct=ci;
      ct.refine(v, tries);
      double temp=ct.llh();
      cout<<"-\t"<<temp<<" "<<ct<<endl;
      if(temp<best) best=temp, cf=ct;
    }

    ci=cf; bag.insert(make_pair(best, ci));
    cout<<k<<"\t"<<best<<" "<<ci<<endl;
  }
  ci=bag.begin()->second;

  if(flasher){
    int num=0;
    double r2=0, r[3]={0};
    for(multimap<double, cascade>::const_iterator i=bag.begin(); num<upto && i!=bag.end(); ++i, num++){
      const cascade & c = i->second;
      r[0]+=c.x, r[1]+=c.y, r[2]+=c.z;
      r2+=c.x*c.x+c.y*c.y+c.z*c.z;
    }

    {
      double r2a=0;
      for(int i=0; i<3; i++){
	r[i]/=num;
	r2a+=r[i]*r[i];
      }

      dR=sqrt(r2/num-r2a);
    }
  }
  else if(!mlpd){
    int num=0;
    double rn=0, r2=0, r[3]={0}, n[3]={0};
    for(multimap<double, cascade>::const_iterator i=bag.begin(); num<upto && i!=bag.end(); ++i, num++){
      const cascade & c = i->second;
      for(int i=0; i<3; i++) n[i]+=c.n[i];
      r[0]+=c.x, r[1]+=c.y, r[2]+=c.z;
      r2+=c.x*c.x+c.y*c.y+c.z*c.z;
      rn+=c.x*c.n[0]+c.y*c.n[1]+c.z*c.n[2];
    }

    {
      double rna=0, r2a=0, n2a=0;
      for(int i=0; i<3; i++){
	r[i]/=num, n[i]/=num;
	rna+=r[i]*n[i], r2a+=r[i]*r[i], n2a+=n[i]*n[i];
      }
      rn/=num; rr=n2a<1?(rna-rn)/(1-n2a):0;

      dR=sqrt(r2/num-r2a+2*rr*(rn-rna)+rr*rr*(1-n2a));
      dA=n2a<1?1/(sqrt(sqrt(n2a)*(3-n2a)/(1-n2a))*cv):0;
    }
  }
  else{
    int num=0;
    double a[3][3]={0}, b[3]={0};

    for(multimap<double, cascade>::const_iterator i=bag.begin(); num<upto && i!=bag.end(); ++i, num++){
      const cascade & c = i->second;
      double dot=c.x*c.n[0]+c.y*c.n[1]+c.z*c.n[2];

      a[0][0]+=1-c.n[0]*c.n[0]; a[0][1]-=c.n[1]*c.n[0]; a[0][2]-=c.n[2]*c.n[0]; b[0]+=c.x-dot*c.n[0];
      a[1][0]-=c.n[0]*c.n[1]; a[1][1]+=1-c.n[1]*c.n[1]; a[1][2]-=c.n[2]*c.n[1]; b[1]+=c.y-dot*c.n[1];
      a[2][0]-=c.n[0]*c.n[2]; a[2][1]-=c.n[1]*c.n[2]; a[2][2]+=1-c.n[2]*c.n[2]; b[2]+=c.z-dot*c.n[2];
    }

    gsl_matrix * A = gsl_matrix_alloc(3, 3);
    gsl_vector * B = gsl_vector_alloc(3);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++) gsl_matrix_set(A, i, j, a[i][j]);
      gsl_vector_set(B, i, b[i]);
    }

    gsl_vector * X = gsl_vector_alloc(3);

    gsl_linalg_HH_solve(A, B, X);

    gsl_matrix_free(A);
    gsl_vector_free(B);

    double r[3];
    for(int i=0; i<3; i++) r[i]=gsl_vector_get(X, i);
    gsl_vector_free(X);

    const bool verbose=false;
    if(verbose){ cout<<"X:"; for(int i=0; i<3; i++) cout<<"\t"<<r[i]; cout<<endl; }

    num=0;
    double r2a=0, n[3]={0};
    for(multimap<double, cascade>::const_iterator i=bag.begin(); num<upto && i!=bag.end(); ++i, num++){
      const cascade & c = i->second;
      for(int i=0; i<3; i++) n[i]+=c.n[i];

      double dx=r[0]-c.x, dy=r[1]-c.y, dz=r[2]-c.z;
      double dot=dx*c.n[0]+dy*c.n[1]+dz*c.n[2];
      dx-=dot*c.n[0], dy-=dot*c.n[1], dz-=dot*c.n[2];

      r2a+=dx*dx+dy*dy+dz*dz;
    }

    {
      double n2a=0;
      for(int i=0; i<3; i++){
	n[i]/=num; n2a+=n[i]*n[i];
      }

      dR=sqrt(r2a/num); rr=0;
      dA=n2a<1?1/(sqrt(sqrt(n2a)*(3-n2a)/(1-n2a))*cv):0;
      // for(int i=0; i<3; i++) COG[i]=r[i];
    }
  }

  step();
}

void spsa(cascade & ci){  // simultaneous perturbation stochastic approximation
  const int steps=loop, start=1;
  const double dL0=40;

  for(int k=0; k<steps; k++){
    double an=start/double(start+k), cn=pow(an, 1./4);
    double dr=cn*dR, da=cn*dA, di=cn*dI;

    ci.setp(); ci.setu();
    cascade c1=ci, c2=ci;

    double nz=2*xrnd()-1;
    double nr=sqrt(1-nz*nz), np=2*M_PI*xrnd();
    double nx=nr*cos(np), ny=nr*sin(np);

    c1.x-=dr*nx, c1.y-=dr*ny, c1.z-=dr*nz; c1.rotate(-cv*da);
    c2.x+=dr*nx, c2.y+=dr*ny, c2.z+=dr*nz; c2.rotate( cv*da);

    double as=2*M_PI*xrnd();
    double is=cos(as), ia=sin(as);
    c1.s/=exp(di*is), c1.a/=exp(di*ia);
    c2.s*=exp(di*is), c2.a*=exp(di*ia);

    double l1=c1.llh(), l2=c2.llh();
    double dL=l2-l1, llh=(l1+l2)/2;

    double fr=an*dL/(2*dL0*cn); if(fabs(fr)>cn) fr*=cn/fabs(fr);

    dr=fr*dR, da=fr*dA, di=fr*dI;
    ci.x-=dr*nx, ci.y-=dr*ny, ci.z-=dr*nz; ci.rotate(-cv*da);
    ci.s/=exp(di*is), ci.a/=exp(di*ia);

    llh=ci.llh();
    cout<<k<<"\t"<<llh<<" "<<ci<<endl;
  }
}

double pv(double y1, double y2, double y3){
  double x1=-1, x2=1, x3=0;
  double D=(x1-x2)*(x1-x3)*(x2-x3);
  double A=(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/D;

  if(A<=0) return y1<y2?-1:y2<y1?1:0;
  else{
    double B=(x3*x3*(y1-y2)+x2*x2*(y3-y1)+x1*x1*(y2-y3))/D;
    //  C=(x2*x3*(x2-x3)*y1+x3*x1*(x3-x1)*y2+x1*x2*(x1-x2)*y3)/D;
    double x=-B/(2*A);
    return x<-1?-1:x>1?1:x;
  }
}

void spsawh(cascade & ci){  // simultaneous perturbation stochastic approximation with Hessian
  const int steps=loop, start=1;
  double llh=ci.llh();

  for(int k=0; k<steps; k++){
    double an=start/double(start+k), cn=pow(an, 0.25);

    ci.setp();
    cascade c1=ci, c2=ci;

    c1.move(cn*dR, cn*dA, cn*dI);
    c2.x=2*ci.x-c1.x, c2.y=2*ci.y-c1.y, c2.z=2*ci.z-c1.z;
    double dot=ci.n[0]*c1.n[0]+ci.n[1]*c1.n[1]+ci.n[2]*c1.n[2];
    for(int i=0; i<3; i++) c2.n[i]=2*dot*ci.n[i]-c1.n[i];
    c2.s=ci.s*ci.s/c1.s; c2.a=ci.a*ci.a/c1.a;

    double l1=c1.llh(), l2=c2.llh();

    double fr=pv(l1, l2, llh)*(an/cn);
    // cout<<l1<<" "<<llh<<" "<<l2<<" "<<fr<<endl<<c1<<endl<<ci<<endl<<c2<<endl;

    if(fr<0){
      ci.x=-fr*c1.x+(1+fr)*ci.x; ci.y=-fr*c1.y+(1+fr)*ci.y; ci.z=-fr*c1.z+(1+fr)*ci.z;
      for(int i=0; i<3; i++) ci.n[i]=-fr*c1.n[i]+(1+fr)*ci.n[i];
      ci.s=pow(c1.s, -fr)*pow(ci.s, 1+fr); ci.a=pow(c1.a, -fr)*pow(ci.a, 1+fr);
    }
    else if(fr>0){
      ci.x=fr*c2.x+(1-fr)*ci.x; ci.y=fr*c2.y+(1-fr)*ci.y; ci.z=fr*c2.z+(1-fr)*ci.z;
      for(int i=0; i<3; i++) ci.n[i]=fr*c2.n[i]+(1-fr)*ci.n[i];
      ci.s=pow(c2.s, fr)*pow(ci.s, 1-fr); ci.a=pow(c2.a, fr)*pow(ci.a, 1-fr);
    }
    ci.normalize();

    llh=ci.llh();
    cout<<k<<"\t"<<llh<<" "<<ci<<endl;
  }
}

void half(){
  dR*=M_SQRT1_2, dA*=M_SQRT1_2;
}

void mcmc(cascade & ci){  // Markov chain monte carlo
  const int steps=loop, upto=steps/2;
  for(int k=0, n=0; k<steps; k++){
    double best=ci.llh();
    ci.setp();

    cascade cf=ci; cf.move();
    double temp=cf.llh();

    if(temp<best) ci=cf, best=temp, n=0;
    else if(k<upto){
      n++;
      if(n>25) half(), step(), n=0;
    }

    cout<<k<<"\t"<<best<<" "<<ci<<endl;
  }
}

void errors(cascade & ci){  // uniform sample to estimate uncertainties
  const int steps=loop;
  ci.setp();

  for(int k=0, n=0; k<steps; k++){
    double best=ci.llh();
    cout<<k<<"\t"<<best<<" "<<ci<<endl;

    cascade cf=ci; cf.move();
    double temp=cf.llh();

    cout<<k<<"\t"<<temp<<" "<<cf<<endl;
  }
}

void err_mc(cascade & ci){  // uniform sample to estimate uncertainties
  const int steps=10*loop;
  ci.setp();

  for(int k=0, n=0; k<steps; k++){
    cascade cf=ci; cf.move();
    double aux=cf.llh();

    if(aux<LX){
      ci=cf; ci.setp();
      cout<<"+";
    }
    else{
      cout<<"-";
    }

    cout<<"\t"<<aux<<" "<<cf<<endl;
  }
}

void method(cascade & ci, int n){
  cout<<"Method: "<<n<<endl;
  if(n<0){
    double llh=0;
    vector< pair<xkey, hix> > sav[cnum];

    switch(n){
    case -1: llh=ci.llh(sav, true,  true); break;
    case -2: llh=ci.llh(sav, false, true); break;
    case -3: llh=ci.llh(sav, true, false); break;
    }

    cout<<"*\t"<<llh<<" "<<ci<<endl;
    ci.printout(sav, cnum);
  }
  else if(n<10){
    double llh=0;
    switch(n){
    case 0: llh=ci.llh(); break;
    case 1: llh=ci.llh(NULL, true,  true); break;
    case 2: llh=ci.llh(NULL, false, true); break;
    case 3: llh=ci.llh(NULL, true, false); break;
    }

    cout<<"?\t"<<llh<<" "<<ci<<endl;
  }
  else{
    switch(n){
    case 10: method(ci, 0); method(ci, 11); break;
    case 11: for(int i=0; i<10; i++) loc_rs(ci); break;
    case 12: spsa(ci);   break;
    case 13: spsawh(ci); break;
    case 14: mcmc(ci);   break;
    case 15: errors(ci); break;
    case 16: err_mc(ci); break;
    }
  }
}

int main(int arg_c, char *arg_a[]){
  // srand48(time(0));
  gsl_set_error_handler_off();

  cascade ci;

  {
    char * bmp=getenv("FSEP");
    if(bmp!=NULL) fsep=atoi(bmp);
    cerr<<"FSEP="<<fsep<<endl;
  }
  {
    char * bmp=getenv("SREP");
    if(bmp!=NULL) srep=atoi(bmp);
    cerr<<"SREP="<<srep<<endl;
  }
  {
    char * bmp=getenv("DREP");
    if(bmp!=NULL) drep=atoi(bmp);
    cerr<<"DREP="<<drep<<endl;
  }
  {
    char * bmp=getenv("LOOP");
    if(bmp!=NULL) loop=atoi(bmp);
    cerr<<"LOOP="<<loop<<endl;
  }
  {
    char * bmp=getenv("NORM");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)) norm=true;
    cerr<<"Normalize everything: "<<(norm?"enabled":"not set")<<"."<<endl;
  }
  {
    char * bmp=getenv("SUMW");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)) sumw=true;
    cerr<<"Multi-PMT sensors: "<<(sumw?"combine PMTs":"treat PMTs individually")<<endl;
  }
  {
    char * bmp=getenv("FAIL");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)) fail=true;
    cerr<<"Fail on warnings: "<<(fail?"enabled":"not set")<<"."<<endl;
  }
  {
    char * bmp=getenv("FAST");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)) fast=true;
    cerr<<"Running in a "<<(fast?"fast":"slow")<<" mode!"<<endl;
  }
  {
    char * bmp=getenv("MLPD");
    if(bmp!=NULL){
      if(*bmp==0 || atoi(bmp)<0) oldf=true, munf=0;
      else oldf=false, munf=atoi(bmp); mlpd=munf>0;
    }
    cerr<<"Running in a "<<(oldf?"compatibility":mlpd?"millipede":"cascade")<<" mode!"<<endl;
  }
  {
    char * bmp=getenv("FLSH");
    if(bmp!=NULL){
      fla.first=strtol(bmp, &bmp, 10);
      fla.second=strtol(bmp+1, &bmp, 10);
      flasher=true;
      cerr<<"Flasher at "<<fla.first<<","<<fla.second<<"."<<endl;
    }
  }
  {
    char * bmp=getenv("FDUR");
    if(bmp!=NULL) fdur=atof(bmp);
    cerr<<"FDUR="<<fdur<<endl;
  }
  {
    char * bmp=getenv("QSAT");
    if(bmp!=NULL) qsat=atof(bmp);
    cerr<<"QSAT="<<qsat<<endl;
  }
  {
    int tmp=0;
    char * bmp=getenv("CNUM");
    if(bmp!=NULL) tmp=atoi(bmp);
    if(tmp>0) cnum=tmp;
    cerr<<"CNUM="<<cnum<<endl;

    cbin=2./cnum;
    xppc::set_res(0.1f, cbin);
  }
  {
    char * bmp=getenv("LSIG");
    if(bmp!=NULL) LSSL::SIG=atof(bmp);
    cerr<<"LSIG="<<LSSL::SIG<<endl;
  }

  xppc::start();
  xppc::initialize(1.f); // flasher?1:0.9*0.85
  xppc::choose(-1);
  xppc::ini();

  bool flag=ci.ini();
  if(flasher){
    xppc::DOM f=xppc::flset(fla.first, fla.second);
    if(!flag) ci.x=f.r[0], ci.y=f.r[1], ci.z=f.r[2];
  }

  cerr<<"Initialized!"<<endl;

  if(flasher && mlpd){
    if(munf%10==3){
      const float * dir = xppc::fldir();
      for(int i=0; i<3; i++) fln[i]=dir[i];
      cerr<<"LED direction at "<<fln[0]<<" "<<fln[1]<<" "<<fln[2]<<endl;
    }

    char * bmp=getenv("FLOR");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)){
      xppc::ikey om; om.str=abs(fla.first), om.dom=fla.second;
      if(ang.cx.find(om)!=ang.cx.end()){
	V<3>& n = ang.cx[om];
	for(int i=0; i<3; i++) fln[i]=n[i];
	cerr<<"Flasher board will be tilted with n=("<<n[0]<<","<<n[1]<<","<<n[2]<<")"<<endl;
      }
    }
  }

  if(arg_c>1){
    for(int i=1; i<arg_c; i++) method(ci, atoi(arg_a[i]));
  }
  else if(arg_c>0){
    cerr<<"Use as following:"<<endl;
    cerr<<"[FLSH=dom,str] [MLPD=0/1] [SREP=sim. events] [DREP=dat. events] "<<arg_a[0]<<" [method]"<<endl;
  }

  if(mlpd || !flasher) xppc::fin();
  xppc::stop();
}
