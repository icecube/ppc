#define LMAX 80      // number of dust loggers
#define LYRS 172     // number of depth points

#define CTX  11
#define CTY  10
#define DIR1 9.3f
#define DIR2 129.3f

#define CX   21
#define CY   19
#define NSTR 240

#define OVER 128     // minimum number of photons per particle (optimized)
#define NTOT 1024
#define NPHO 1024  // maximum number of photons propagated by one thread

#define WNUM   512   // number of wavelength slices
#define MAXLYS 172   // maximum number of ice layers
#define MAXGEO 16384 // maximum number of OMs
#define MAXRND 131072 // max. number of random number multipliers

// #define DTMN
#define XXX 1.e-5f
#define FPI 3.141592653589793f
#define OMR 0.16510f  // DOM radius [m]

static const float doma=FPI*OMR*OMR; // DOM cross-sectional area [m^2]
static const float omav=0.335822;    // average DOM angular (lab) sensitivity
static const float dppm=2450.08;     // photons per meter for nominal IceCube DOM
static const float fcv=FPI/180.f;

static unsigned int ovr=1;

float xrnd();

struct DOM{
  float R, F;
  float r[3];
};

struct ikey{
  int str, dom;

  bool isinice() const{
    return str>86 || (str>0 && dom>=1 && dom<=60);
  }

  bool operator< (const ikey & rhs) const {
    return str == rhs.str ? dom < rhs.dom : str < rhs.str;
  }

  bool operator!= (const ikey & rhs) const {
    return str != rhs.str || dom != rhs.dom;
  }
};

struct OM:DOM,ikey{};
vector<OM> i3oms;
map<int, pair<float, float> > strs;

template <int n> class V{
  float x[n];

public:
  V(){
    for(int i=0; i<n; i++) x[i]=i<n-1?0:1;
  }

  float & operator[](int i){
    return x[i];
  }

  float dot(V<n> q){
    float sum=0;
    for(int i=0; i<n; i++) sum+=x[i]*q[i];
    return sum;
  }
};

struct name:ikey{
  int omt;  // module type
  int type; // wavelength behavior
  float rde, hv;
  float azi; // cable
  V<3> tilt;

  name(){}
  name(ikey k, int m, int t, float r, float h, V<3>& axis, float ph):ikey(k){
    omt=m; type=t; rde=r; hv=h;
    tilt=axis; azi=ph;
  }
};

struct mesh{
  vector< V<3> > dirs;

  int ini(const string & file){
    ifstream inFile(file.c_str(), ifstream::in);
    if(!inFile.fail()){
      string in;
      while(getline(inFile, in)){
	V<3> dir;
	if(3==sscanf(in.c_str(), "%*d %f %f %f", &dir[0], &dir[1], &dir[2])){
	  dirs.push_back(dir);
	}
      }
      inFile.close();
      cerr<<"Configured "<<dirs.size()<<" uniform directions"<<endl;
    }
    return dirs.size();
  }
} ico;

bool nextgen=false;

struct itype{
  float area, beta, rde, fx, Rr, Rz, cable;
  vector< V<3> > dirs;

  bool def;
  float ave;  // average angular sensitivity
  float mas;  // maximum angular sensitivity
  vector<float> s; // ang. sens. coefficients

  itype(): def(false){ }

  void add(string file){
    def=true; mas=1, ave=0;
    area=1, beta=0.33f, Rr=OMR, Rz=OMR, cable=0;

    ifstream inFile(file.c_str(), ifstream::in);
    if(!inFile.fail()){
      float aux;
      if(inFile >> aux) mas=aux;
      while(inFile >> aux) s.push_back(aux);
      if(s.empty()) s.push_back(1.f);

      if(mas>0){
	ave=s[0];
	for(unsigned int i=2; i<s.size(); i+=2) ave+=s[i]/i;
      }
      else{
	ave=(1-s[0])/2;
      }
      if(ave>0) cerr<<"Loaded "<<(s.size()+1)<<" angsens coefficients"<<endl;
      else{ cerr<<"File "<<file<<" did not contain valid data"<<endl; exit(1); }
      inFile.close();
    }
    else{ cerr<<"Could not open file "<<file<<endl; exit(1); }

    ave=omav; // for backwards compatibility
    rde=(mas>0?mas:1.f)/(doma*ave);
  }

  void add(float th, float ph){
    float ct=cos(th*fcv), st=sin(th*fcv);
    float cp=cos(ph*fcv), sp=sin(ph*fcv);
    V<3> dir;
    dir[0]=st*cp, dir[1]=st*sp, dir[2]=ct;
    dirs.push_back(dir);
  }

  float aS(float x){
    float al=acos(x);
    return al-sin(2*al)/2;
  }

  float f(float x){ // angular sensitivity curve (peaks at 1)
    if(beta<-1) return sqrt(1-x*x);
    else{
      float sum=x>0?x:0;
      if(beta<1){
	float y=sqrt(1-x*x)/beta;
	if(y>1){
	  y=1/y;
	  float c=1-beta*beta;
	  sum+=(aS(y)/c-aS(y*fabs(x)/sqrt(c))*fabs(x))/FPI;
	}
      }
      return sum;
    }
  }

  float xarea(float dot){ // OM cross-sectional area for direction with cos(zenith)=dot
    return Rz>0?FPI*Rr*sqrt(Rz*Rz-dot*dot*(Rz*Rz-Rr*Rr)):-4*Rr*Rz*sqrt(1-dot*dot);
  }

  void fraq(){ // max (over all directions) of PMT to sensor cross-sectional area
    if(def) return;

    float fr_flat=0, sum_ave=0;
    for(vector< V<3> >::iterator i=ico.dirs.begin(); i!=ico.dirs.end(); ++i){
      float tot=0;
      for(vector< V<3> >::const_iterator j=dirs.begin(); j!=dirs.end(); ++j) tot+=f(i->dot(*j));
      float fra=xarea((*i)[2]); fra=fra>0?tot/fra:0;
      if(fr_flat<fra) fr_flat=fra; sum_ave+=tot;
    }
    sum_ave/=ico.dirs.size();

    fx=fr_flat;
    rde=area*fr_flat/sum_ave;
  }

  int getPMT(V<3> dir, V<3> pos, V<3> tilt, float rnd, float & ht, float ph = -1.f){
    if(beta==-2){ // model from code and thesis of Jakob Beise
      const float cn=0.299792458/1.46; // propagation speed in quartz
      const float abs=3.0/cn; // absorption length in ns
      const float att=0.678; // attenuation length along cylinder
      const float cut=0.01; // cut in sigma gamma parameterization
      const float tau=1.6; // paint re-emission decay time
      const float h=0.76; // inner tube height in m
      const float a=h/(1-exp(-h/att)); // effective att. at distance h

      float z=Rr*pos[2]; // z-coordinate from midpoint
      z*=h/fabs(2*Rz); // parameterizaion only for WOM inner tube height of 76 cm

      int i=rnd-0.5<z/(2*a)?1:0;
      z=i>0?h/2-z:h/2+z; // distance to detection PMT

      float d=(z/h)*(1-z/(2*a))/(1-h/(2*a)); // simulating acceptance z-dependence by
      z=-att*log(1-d*(1-exp(-h/att))); // shifting z slightly (up to -1.2 ... +0.4 cm)

      if(z<cut) z=cut; // paramerized only down to 1 cm
      float x=log(z), dt;
      float s=exp(((0.0034112*x+0.041727)*x+0.45550)*x-1.9188); // sigma
      float g=exp(((0.0012658*x+0.017228)*x-0.012750)*x+0.20924); // gamma
      do{
	dt=s*sqrt(2*g*pow(1-xrnd(), 1/(1-g))-1); // King function (transit pdf)
      } while(dt>-abs*log(xrnd())); // rejecting absorbed photons, otherwise tail is too long
      dt+=z/cn-tau*log(xrnd()); // adding tz (shortest distance) and paint re-emission time

      ht+=dt;
      return i;
    }
    else if(def){
      bool flag;
      if(mas>0){
	float sum;
	{
	  float x = dir.dot(tilt);
	  float y=1;
	  sum=s[0];
	  for(unsigned int i=1; i<s.size(); i++){ y*=x; sum+=s[i]*y; }
	}

	flag=mas*rnd<sum;
      }
      else{
	flag=pos.dot(tilt)>s[0];
      }
      return flag?0:-1;
    }
    else{
      if(ph>=0){ // rotating photon by -ph instead of PMTs
	ph-=cable;
	float cp=cos(fcv*ph);
	float sp=sin(fcv*ph);
	float nx=dir[0]*cp+dir[1]*sp;
	float ny=dir[1]*cp-dir[0]*sp;
	dir[0]=nx, dir[1]=ny;
      }

      float crx=fx*xarea(dir[2]);

      unsigned int k=0;
      float tot=0;
      for(vector< V<3> >::const_iterator j=dirs.begin(); j!=dirs.end(); ++j, k++){
	tot+=f(-dir.dot(*j))/crx;
	if(tot>rnd) break;
      }

      return k<dirs.size()?(int) k:-1;
    }
  }
};

struct spec{
  bool swv; // single wavelength distribution
  int num;  // number of kept wavelength bins
  float bin, wmin, wmax;

  spec(): swv(false), num(450), wmin(250.f), wmax(700.f){
    bin=(wmax-wmin)/num;
  }

  void ssw(float wva){ // set single wavelength
    swv=true; num=10; wmin=wva-5.f; wmax=wva+5.f;
    bin=(wmax-wmin)/num;
  }

  float wav(float i){
    return wmin+bin*i;
  }

  float np(float wv, float * ng = NULL){ // phase and group refrative indices
    float np=1.55749-wv*(1.57988-wv*(3.99993-wv*(4.68271-2.09354*wv)));
    if(ng!=NULL) *ng=np*(1+0.227106-wv*(0.954648-wv*(1.42568-0.711832*wv)));
    return np;
  }

  float cherenkov(float wva){
    if(swv){ // over 10 nm
      return 1/(wmax-wmin);
    }
    else{ // um^-2 (here) * nm*cm^2 (later) = 0.1 m
      const float c=FPI*0.2/137.03599911;
      float wv=wva*1.e-3, n=np(wv);
      return c*(1-1/(n*n))/(wv*wv);
    }
  }
} qwv;

struct irde{
  float rde;
  vector<float> dat, sum, rat;

  int oms;  // number of OMs of this type
  float rmax; // max RDE for this type

  irde(): rde(0.f), oms(0), rmax(0.f){}

  void addr(float r){ // enter new OM
    if(r>rmax) rmax=r;
    oms++;
  }

  float binf(float p){
    int i, j=0; while(sum[++j]<=0); j--; while(sum[++j]<p); i=j-1;
    return i + (sum[j]>sum[i] ? (p-sum[i])/(sum[j]-sum[i]) : 0);
  }

  void read(const vector<float> & qw, const vector<float> & qf, float blo = -1.f, float bhi = -1.f){
    dat.clear();

    int k=0, n=qw.size();
    float wlo=qw[0]-(blo<0?qw[1]-qw[0]:blo);
    float whi=qw[n-1]+(bhi<0?qw[n-1]-qw[n-2]:bhi);
    float w1=wlo, w2=qw[0], f1=0, f2=qf[0];

    for(float w=qwv.wmin+qwv.bin/2; w<qwv.wmax; w+=qwv.bin){
      while(qw[k]<w && k<n){
	w1=w2, f1=f2; k++;
	if(k<n) w2=qw[k], f2=qf[k];
	else w2=whi, f2=0;
      }
      dat.push_back(wlo<w && w<whi ? (f1*(w2-w)+f2*(w-w1))/(w2-w1) : 0);
    }
  }
};

map<ikey, int> omts;
map<int, itype> types;
map<pair<int,int>, irde> irdes;
irde env;

map<ikey, float> hvs;
map<ikey, pair<float, int> > rdes;
map<ikey, V<3> > cx;
map<ikey, float > dx;

struct hit{
  unsigned int i, n, z;
  float t;
  float pth, pph, dth, dph;
};

#ifdef XCPU
struct int4{
  int x, y, z, w;
};

struct float2{
  float x, y;
};

struct float3:float2{
  float z;
};

struct float4:float3{
  float w;
};
#endif

struct pbuf{
  float4 r;        // location, time
  float4 n;        // direction
  unsigned int q;  // wavelength bin
  unsigned int i;  // segment number
  int fla, ofla;
};

struct photon{
  float4 r;    // location, time
  float4 n;    // direction, track length
  unsigned int q; // track segment
  unsigned int num; // number of photons in this bunch
  int type;    // source type
  float f;     // fraction of light from muon alone (without cascades)
  union{
    struct{
      float a, b;  // longitudinal development parametrization coefficients
      float beta;  // velocity of incident particle
      float tau;   // luminescence decay time
    };
    struct{
      float ka, up; // 2d-gaussian rms and zenith of cone
      float fldr;   // horizontal direction of the flasher led #1
      short fla, ofla;
    };
    int4 c;  // for quick copy
  };
};

struct ices{
  int wvl;               // wavelength number of this block
  float ocm;             // 1 / speed of light in medium
  float coschr, sinchr;  // cos and sin of the cherenkov angle
  struct{
    float abs;           // absorption
    float sca;           // scattering
  } z [MAXLYS];
};

struct aniz{
  float k1;
  float k2;
  float ra;
  float rb;
};

struct line{
  short n, max;
  float x, y, r;
  float h, d;
  float dl, dh;
};

struct datz{
  ices w[WNUM];
  float2 lp[LMAX][LYRS];
  unsigned int rm[MAXRND];
  unsigned long long rs[MAXRND];
};

struct dats{
  atomic<unsigned int> hidx;
  unsigned int ntot;

  float rx;
  float hifl;

  unsigned int hnum;    // size of hits buffer
  int size;   // size of kurt table
  int rsize;  // count of multipliers
  int gsize;  // count of initialized OMs

  short tmod; // tilt model: 1: 1d, 2: 2d
  short vthk; // 0: uniform layers, 1: variable thickness
  float dh, rdh, hmin; // step, 1/step, and min depth

  float ocv;  // 1 / speed of light in vacuum
  float sf;   // scattering function: 0=HG; 1=SAM
  float g, gr; // g=<cos(scattering angle)> and gr=(1-g)/(1+g)

  float xR;   // DOM oversize scaling factor
  float SF, G, GR; // hole ice sf, g, gr
  float hr, hr2, hs, ha; // hole ice radius, radius^2, effective scattering and absorption coefficients

  float azx, azy;  // ice anisotropy direction

  int cn[2];
  float cl[2], crst[2];

  float cb[2][2];

  int lnum, lpts, l0;
  float lmin, lrdz, r0;
  float lnx, lny;
  float lr[LMAX];

  float mmin[2], mstp[2];
  int mnum[2], mcut[2];

  float k1, k2, kz, fr; // ice absorption anisotropy parameters
  aniz az[MAXLYS];

  unsigned short ls[NSTR];
  unsigned short is[CX][CY];
  char mcol[CTX][CTY];

  float sum, bfr[12];

  line sc[NSTR];
};

struct doms{
  DOM oms[MAXGEO];
  name names[MAXGEO];
  struct{
    float w, i, f;

    float x(){
      return i+(f-i)*xrnd();
    }
  } wvs [WNUM];

  hit * hits;
  photon * pz;
  pbuf * bf;

  float eff;  // OM efficiency correction
} q;

dats * ed;
datz * ez;
DOM * oms;

unsigned short sname(int n){
  static map<int, unsigned short> overflow;
  name & s = q.names[n];
  if(s.str>86){
    static unsigned short next_num=86+10+1;
    map<int, unsigned short>::iterator it=overflow.find(s.str);
    if(it==overflow.end()){
      if(next_num>=0x8000){ cerr<<"Number of string labels exceeds capacity of "<<(unsigned int) 0x8000<<endl; exit(1); }
      overflow.insert(make_pair(s.str, next_num++));
    }
    return overflow[s.str];
  }
  else return (unsigned short)( (s.str>78 && s.dom>10) ? s.str+10 : s.str );
}

static const float zoff=1948.07;
unsigned int sv=0;

void rs_ini(){
  dats & d = * ed;
  datz & z = * ez;

  union{
    unsigned long long da;
    struct{
      unsigned int ax;
      unsigned int dx;
    };
  } s;

  s.ax=362436069;
  s.dx=1234567;

  s.ax+=sv;

  for(int i=0; i<d.rsize; i++) z.rs[i]=s.da;
}

struct ini{
  float ctr(line & s, int m){
    dats & d = * ed;
    return d.cb[0][m]*s.x+d.cb[1][m]*s.y;
  }

  ini(){
    ed = new dats;
    ez = new datz;
  }

  ~ini(){
    delete ed;
    delete ez;
  }

  void set(){
    dats & d = * ed;
    datz & z = * ez;

    {
      d.hidx=0;
    }

    string ppcdir("");
    {
      char * env = getenv("PPCTABLESDIR");
      if(env!=NULL) ppcdir=string(env)+"/";
      else{
	env = getenv("I3_SRC");
	if(env!=NULL) ppcdir=string(env)+"/ppc/resources/ice/";
      }
      cerr<<"Configuring ppc in \""<<ppcdir<<"\""<<endl;
    }

    string omdir("");
    {
      char * env = getenv("NEXTGENDIR");
      if(env!=NULL) omdir=string(env)+"/";
      else omdir=ppcdir;
      cerr<<"Configuring nextgen sensors in \""<<omdir<<"\""<<endl;
    }

    string icedir("");
    {
      char * env = getenv("ICEMODELDIR");
      if(env!=NULL) icedir=string(env)+"/";
      else icedir=ppcdir;
      cerr<<"Configuring icemodel in \""<<icedir<<"\""<<endl;
    }

    string tiltdir("");
    {
      char * env = getenv("TILTMODELDIR");
      if(env!=NULL) tiltdir=string(env)+"/";
      else tiltdir=icedir;
      cerr<<"Configuring tiltmodel in \""<<tiltdir<<"\""<<endl;
    }

    string holeice("");
    {
      char * env = getenv("PPCHOLEICE");
      if(env!=NULL) holeice=string(env);
      else holeice=ppcdir+"as.dat";
      cerr<<"Configuring holeice from \""<<holeice<<"\""<<endl;
    }

    float dk1, dk2, dkz; // ice anisotropy parameters

    {
      ifstream inFile((icedir+"cfg.txt").c_str(), ifstream::in);
      if(!inFile.fail()){
	string in;
	float aux;
	vector<float> v;
	while(getline(inFile, in)) if(sscanf(in.c_str(), "%f", &aux)==1) v.push_back(aux);

	{
	  char * OVSZ=getenv("OVSZ");

	  if(OVSZ!=NULL){
	    float ovsz=atof(OVSZ); if(v.size()>=1) v[0]=ovsz;
	    cerr<<"Using oversize as set with OVSZ="<<ovsz<<endl;
	  }
	}

	if(v.size()>=4){
	  int xR=lroundf(v[0]); d.xR=xR; ovr*=xR*xR;
	  q.eff=v[1], d.sf=v[2], d.g=v[3]; d.gr=(1-d.g)/(1+d.g);
	  cerr<<"Configured: xR="<<xR<<" eff="<<q.eff<<" sf="<<d.sf<<" g="<<d.g<<endl;

	  if(v.size()<12) d.SF=d.sf, d.G=d.g, d.GR=d.gr;
	  else d.SF=v[10], d.G=v[11], d.GR=(1-d.G)/(1+d.G);

	  if(v.size()>=10){
	    float xH=v[7], hS=v[8], hA=v[9];
	    d.hr=OMR*xH; d.hr2=d.hr*d.hr; d.hs=1/(hS*(1-d.G)), d.ha=1/hA;
	    if(xH>0) cerr<<"With hole ice: xH="<<xH<<" sca="<<hS<<" ("<<d.SF<<","<<d.G<<") abs="<<hA<<endl;
	  }
	  else d.hr=0, d.hr2=0, d.hs=0, d.ha=0;

	  if(v.size()>=7){
	    const float thx=v[4];
	    d.azx=cos(fcv*thx), d.azy=sin(fcv*thx);
	    dk1=exp(v[5]); dk2=exp(v[6]); dkz=1/(dk1*dk2);
	    cerr<<"Ice anisotropy is k("<<thx<<")="<<dk1<<","<<dk2<<","<<dkz<<endl;
	  }
	  else dk1=1, dk2=1, dkz=1, d.azx=1, d.azy=0;

	  if(v.size()>=15){
	    // new absorption anisotropy
	    d.k1=exp(v[12]); d.k2=exp(v[13]); d.kz=exp(v[14]);
	    cerr<<"New Ice anisotropy is "<<d.k1<<","<<d.k2<<","<<d.kz<<endl;
	    if(d.k1>=d.k2 && d.k2==d.kz){
	      float r=d.k1/d.k2;
	      float s=sqrt(r*r-1);
	      r=(s>XXX?log(r+s)/s:1+s/2)/d.k2;
	      d.k1*=r, d.k2*=r, d.kz*=r;
	      cerr<<"Renorm. NI anisotropy "<<d.k1<<","<<d.k2<<","<<d.kz<<endl;
	    }
	    else{
	      cerr<<"Warning: taking absorption anisotropy as given (not renormalizing)!"<<endl;
	    }
	  }
	  else d.k1=1, d.k2=1, d.kz=1;

	  if(v.size()>=16){
	    // scaling for absorption anisotropy (old implementation)
	    d.fr=v[15];
	    cerr<<"Ice absorption anisotropy scaling is "<<d.fr<<endl;
	  }
	  else d.fr=1;

	  if(v.size()>=28){
	    for(int i=0; i<12; i++) d.bfr[i]=v[16+i];

	    {
	      char * BFRA=getenv("BFRA");
	      float bfra=BFRA==NULL?1.0:atof(BFRA);

	      char * BFRB=getenv("BFRB");
	      float bfrb=BFRB==NULL?1.0:atof(BFRB);

	      if(BFRA!=NULL || BFRB!=NULL) cerr<<"Setting BFRA="<<bfra<<" BFRB="<<bfrb<<endl;
	      for(int i=0; i<12; i+=4) d.bfr[i]*=i<8?sqrt(bfra):bfra*bfrb;
	    }

	    {
	      float step=0.01, sum=0;
	      for(float x=step/2; x<1; x+=step){
		float y=sqrt(1-x*x);
		float sx=max(0.f, d.bfr[0]*exp(-d.bfr[1]*pow(atan(d.bfr[3]*y), d.bfr[2])));
		float sy=max(0.f, d.bfr[4]*exp(-d.bfr[5]*pow(atan(d.bfr[7]*y), d.bfr[6])));
		float mx=max(0.f, d.bfr[8]*atan(d.bfr[11]*y*x)*exp(-d.bfr[9]*y+d.bfr[10]*x));
		sum+=sx*sx+sy*sy+mx*mx;
	      }
	      sum*=step/2; d.sum=sum;
	    }
	    cerr<<"Initialized BFR diffusion patterns; s_eff="<<d.sum<<" m^-1"<<endl;
	  }
	  else{
	    d.sum=0;
	    for(int i=0; i<12; i++) d.bfr[i]=i%4<1?0:1;
	  }
	}
	else{ cerr<<"File cfg.txt did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else{ cerr<<"Could not open file cfg.txt"<<endl; exit(1); }
    }

    {
      ifstream inFile((omdir+"om.conf").c_str(), ifstream::in);
      if(!inFile.fail()){
	string in;
	while(getline(inFile, in)){
	  int m;
	  unsigned int n;
	  itype t;
	  float th, ph;
	  float other;
	  int read=sscanf(in.c_str(), "%*s %d %f %f %f %f %d %f %f %f", &m, &t.area, &t.beta, &t.Rr, &t.Rz, &n, &th, &ph, &other);
	  t.cable=read>=9?other:0;
	  if(read>=8){
	    t.add(th, ph);
	    for(unsigned int i=1; i<n && getline(inFile, in); i++)
	      if(2==sscanf(in.c_str(), "%f %f %f", &th, &ph, &other)) t.add(th, ph);
	    if(t.dirs.size()==n) types.insert(make_pair(m, t));
	  }
	}
	inFile.close();
      }

      if(!types.empty()){
	if(ico.ini(omdir+"om.dirs")<1){ cerr<<"Error: could not initialize an array of directions"<<endl; exit(1); }
	nextgen=true;
	for(map<int, itype>::iterator j=types.begin(); j!=types.end(); ++j){
	  j->second.fraq();
	  cerr<<" OM Type "<<j->first<<" with "<<j->second.dirs.size()<<" PMTs added ("<<j->second.rde<<")"<<endl;
	}
      }
    }

    if(types.find(-1)==types.end()){ // read angular sensitivity parameters
      types[-1].add(holeice);
    }

    { // initialize random numbers
      int size;
      vector<unsigned int> rx;

      ifstream inFile((ppcdir+"rnd.txt").c_str(), ifstream::in);
      if(!inFile.fail()){
	string in;
	while(getline(inFile, in)){
	  stringstream str(in);
	  unsigned int a;
	  if(str>>a) rx.push_back(a);
	}
	if(rx.size()<1){ cerr<<"File rnd.txt did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
      else{ cerr<<"Could not open file rnd.txt"<<endl; exit(1); }

      size=rx.size();
      if(size>MAXRND){
	cerr<<"Error: too many random multipliers ("<<size<<"), truncating to "<<MAXRND<<endl;
	size=MAXRND;
      }

      cerr<<"Loaded "<<size<<" random multipliers"<<endl;

#ifndef DTMN
      timeval tv; gettimeofday(&tv, NULL);
      sv=1000000*(unsigned long long)tv.tv_sec+tv.tv_usec;
#endif

      d.rsize=size;
      for(int i=0; i<size; i++) z.rm[i]=rx[i];
    }

    {
      ifstream inFile((ppcdir+"geo-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
	if(!i3oms.empty()){
	  i3oms.clear();
	  cerr<<"Warning: overwriting existing geometry!"<<endl;
	}
	OM om;
	string mbid;
	unsigned long long omid;
	while(inFile>>mbid>>hex>>omid>>dec>>om.r[0]>>om.r[1]>>om.r[2]>>om.str>>om.dom){
	  om.R=OMR, om.F=1;
	  om.r[2]+=zoff;
	  i3oms.push_back(om);
	}
	inFile.close();
      }
    }
    {
      char * HIFL=getenv("HIFL");
      float hifl=HIFL==NULL?0.f:atof(HIFL);

      if(HIFL!=NULL) cerr<<"Setting HIFL="<<hifl<<" ns"<<endl;
      d.hifl=hifl;
    }
    {
      map<ikey, pair<float, float> > geco;

      char * bmp=getenv("GECO");
      if(bmp!=NULL){
	ikey om;
	float x, y, xH, hS;
	int num=sscanf(bmp, "%d %d %f %f %f %f", &om.str, &om.dom, &x, &y, &xH, &hS);
	if(num>=4){
	  geco[om]=make_pair(OMR*x, OMR*y);
	  cerr<<"Applying geometry correction to DOM "<<om.str<<","<<om.dom<<" of "<<x<<","<<y<<" x OMR ("<<(100*OMR)<<" cm)"<<endl;
	}
	if(num>=5){
	  d.hr=OMR*xH; d.hr2=d.hr*d.hr;
	  cerr<<"Updating hole ice parameters: xH="<<xH<<endl;
	}
	if(num>=6){
	  d.hs=1/(hS*(1-d.G));
	  cerr<<"Updating hole ice parameters: sca="<<hS<<endl;
	}
      }

      for(vector<OM>::iterator i=i3oms.begin(); i!=i3oms.end(); ++i){
	map<ikey, pair<float, float> >::iterator j=geco.find(*i);
	if(j!=geco.end()){
	  pair<float, float>& gc = j->second;
	  i->r[0]+=gc.first, i->r[1]+=gc.second;
	}
      }
    }
    {
      ifstream inFile((ppcdir+"str-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
	cerr<<"Using fixed positions of hole ice columns from str-f2k!"<<endl;

	int str;
	float x, y;
	while(inFile>>str>>x>>y) strs[str]=make_pair(x, y);
	inFile.close();
      }
    }
    {
      ifstream inFile((ppcdir+"eff-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
	if(!rdes.empty()){
	  rdes.clear();
	  cerr<<"Warning: overwriting existing RDE table!"<<endl;
	}
	int typ;
	ikey om;
	float eff;

	string in;
	while(getline(inFile, in)){
	  int num=sscanf(in.c_str(), "%d %d %f %d", &om.str, &om.dom, &eff, &typ);
	  if(num<4) typ=0;
	  if(num>=3) if(om.isinice()) rdes.insert(make_pair(om, make_pair(eff, typ)));
	}

	inFile.close();
      }
    }
    {
      ifstream inFile((omdir+"om.map").c_str(), ifstream::in);
      if(!inFile.fail()){
	if(!omts.empty()){
	  omts.clear();
	  cerr<<"Warning: overwriting existing OM type table!"<<endl;
	}
	ikey om;
	int m;
	while(inFile>>om.str>>om.dom>>m) if(om.isinice()) omts.insert(make_pair(om, m));
	inFile.close();
      }
    }
    {
      ifstream inFile((ppcdir+"hvs-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
	if(!hvs.empty()){
	  hvs.clear();
	  cerr<<"Warning: overwriting existing HV table!"<<endl;
	}
	ikey om;
	float hv;
	while(inFile>>om.str>>om.dom>>hv) if(om.isinice()) hvs.insert(make_pair(om, hv));
	inFile.close();
      }
    }

    {
      ifstream inFile((ppcdir+"cx.dat").c_str(), ifstream::in);
      if(!inFile.fail()){
	ikey om;
	V<3> dir;
	float r;
	while(inFile >> om.str >> om.dom >> dir[0] >> dir[1] >> dir[2] >> r) cx[om]=dir;
	if(cx.size()>0) cerr<<"Loaded "<<cx.size()<<" DOM orientations"<<endl;
	else{ cerr<<"File cx.dat did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
    }
    {
      ifstream inFile((ppcdir+"dx.dat").c_str(), ifstream::in);
      if(!inFile.fail()){
	ikey om;
	float dir;
	float r;
	while(inFile >> om.str >> om.dom >> dir >> r){
	  while(dir<0) dir+=360.f; // negative value means unset
	  dx[om]=dir;
	}
	if(dx.size()>0) cerr<<"Loaded "<<dx.size()<<" cable positions"<<endl;
	else{ cerr<<"File dx.dat did not contain valid data"<<endl; exit(1); }
	inFile.close();
      }
    }

    float Rr=0, Rz=0;

    { // initialize geometry
      vector<DOM> oms;
      vector<name> names;
      int nhqe=0;

      sort(i3oms.begin(), i3oms.end());
      for(vector<OM>::iterator i=i3oms.begin(); i!=i3oms.end(); ++i) if(i->isinice()){
	ikey om(*i);

	int m, t;
	float r, h;
	V<3> tilt;
	float azi = -1.f;

	if(omts.empty()) m=-1;
	else{
	  map<ikey, int>::iterator j=omts.find(om);
	  m=j==omts.end()?-1:j->second;
	  map<int, itype>::const_iterator it=types.find(m);
	  if(it==types.end()) m=-1;

	  it=types.find(m);
	  if(it!=types.end()){
	    const itype & t=it->second;
	    i->R=t.Rr, i->F=t.Rz/t.Rr;
	  }
	}

	Rr=max(Rr, i->R), Rz=max(Rz, i->R*fabs(i->F));
	oms.push_back(*i);

	{
	  map<ikey, pair<float, int> >::iterator j=rdes.find(om);
	  if(j!=rdes.end()){
	    nhqe++;
	    r=j->second.first;
	    t=j->second.second;
	  }
	  else r=1, t=0;

	  irdes[make_pair(m,t)].addr(r);
	}

	if(hvs.empty()) h=1200;
	else{
	  map<ikey, float>::iterator j=hvs.find(om);
	  h=j==hvs.end()?0:j->second;
	}

	{
	  map<ikey, V<3> >::iterator ci=cx.find(om);
	  if(ci!=cx.end()) tilt=ci->second;

	  map<ikey, float >::iterator di=dx.find(om);
	  if(di!=dx.end()) azi = di->second;
	}

	names.push_back(name(om, m, t, r, h, tilt, azi));
      }
      if(nhqe>0) cerr<<"Loaded "<<nhqe<<" RDE coefficients"<<endl;
      for(map<pair<int,int>, irde>::const_iterator i=irdes.begin(); i!=irdes.end(); ++i)
	cerr<<" Found "<<i->second.oms<<" OMs of type "<<i->first.first<<","<<i->first.second<<endl;

      int gsize = oms.size();
      if(gsize>MAXGEO){
	cerr<<"Error: too many OMs ("<<gsize<<"), truncating to "<<MAXGEO<<endl;
	gsize=MAXGEO;
      }

      for(int n=0; n<gsize; n++){ q.oms[n]=oms[n]; q.names[n]=names[n]; }

      d.gsize=gsize;
    }

    Rr*=d.xR, Rz*=d.xR;

    map<unsigned short, short> num;
    {
      map<unsigned short, float> l;
      map<unsigned short, line> sc;
      for(int n=0; n<d.gsize; n++){
	unsigned short str=sname(n);
	line & s = sc[str];
	DOM & om = q.oms[n];
	if(num.find(str)==num.end()){
	  l[str]=om.r[2];
	  s.h=om.r[2];
	  s.n=n;
	}
	else{
	  if(l[str]>om.r[2]) l[str]=om.r[2];
	  if(s.h<om.r[2]) s.h=om.r[2];
	  if(s.n>n) s.n=n;
	}
	num[str]++; s.x+=om.r[0], s.y+=om.r[1];
      }
      if(sc.size()>NSTR){ cerr<<"Number of strings exceeds capacity of "<<NSTR<<endl; exit(1); }

      for(map<unsigned short, short>::iterator i=num.begin(); i!=num.end(); ++i){
	unsigned short str=i->first, n=i->second;
	line & s = sc[str];
	float d=s.h-l[str];
	if(n>1 && d<=0){ cerr<<"Cannot estimate the spacing along string "<<(int)str<<endl; exit(1); }
	s.x/=n, s.y/=n, s.r=0; s.d=n>1?(n-1)/d:0; s.dl=0; s.dh=0;

	map<int, pair<float, float> >::iterator j=strs.find(str);
	if(j!=strs.end()) s.x=j->second.first, s.y=j->second.second;
      }

      for(int n=0; n<d.gsize; n++){
	unsigned short str=sname(n);
	line & s = sc[str];
	DOM & om = q.oms[n];

	float dx=s.x-om.r[0], dy=s.y-om.r[1];
	float dr=dx*dx+dy*dy; if(dr>s.r) s.r=dr;

	if(s.d>0){
	  float dz=om.r[2]-(s.h+(s.n-n)/s.d);
	  if(s.dl>dz) s.dl=dz; if(s.dh<dz) s.dh=dz;
	}
      }

      d.rx=0;
      int n=0;
      for(map<unsigned short, short>::iterator i=num.begin(); i!=num.end(); ++i, n++){
	unsigned short str=i->first;
	line & s = sc[str];
	s.max=i->second-1;
	i->second=n;
	s.r=Rr+sqrt(s.r);
	if(d.hr>s.r) s.r=d.hr;
	if(d.rx<s.r) d.rx=s.r;
	s.dl-=Rz, s.dh+=Rz;
	d.sc[n]=s;
      }
    }

    float sin12=0;
    {
      float bv[2][2];
      bv[0][0]=cos(fcv*DIR1);
      bv[0][1]=sin(fcv*DIR1);
      bv[1][0]=cos(fcv*DIR2);
      bv[1][1]=sin(fcv*DIR2);

      float det=bv[0][0]*bv[1][1]-bv[0][1]*bv[1][0];
      d.cb[0][0]=bv[1][1]/det;
      d.cb[0][1]=-bv[0][1]/det;
      d.cb[1][0]=-bv[1][0]/det;
      d.cb[1][1]=bv[0][0]/det;

      for(int i=0; i<2; i++) sin12+=bv[0][i]*bv[1][i];

      sin12=sqrt(1-sin12*sin12);
      d.rx/=sin12;
    }

    map<unsigned short, int> cells[CX][CY];
    {
      float cl[2]={0,0}, ch[2]={0,0}, crst[2];

      int n=0;
      for(map<unsigned short, short>::iterator i=num.begin(); i!=num.end(); ++i, n++){
	line & s = d.sc[i->second];
	for(int m=0; m<2; m++){
	  if(n==0 || ctr(s, m)<cl[m]) cl[m]=ctr(s, m);
	  if(n==0 || ctr(s, m)>ch[m]) ch[m]=ctr(s, m);
	}
      }

      d.cn[0]=CX;
      d.cn[1]=CY;

      for(int m=0; m<2; m++){
	float diff=ch[m]-cl[m];
	d.cn[m]=min(d.cn[m], 1+2*(int)lroundf(diff/125));

	if(d.cn[m]<=1){
	  ch[m]=cl[m]=(cl[m]+ch[m])/2;
	  crst[m]=1/(d.rx*(2+XXX)+diff);
	}
	else{
	  float s=Rr*(d.cn[m]-1);
	  if(diff<2*s){
	    cerr<<"Warning: tight string packing in direction "<<(m<1?"x":"y")<<endl;
	    float ave=(cl[m]+ch[m])/2;
	    cl[m]=ave-s; ch[m]=ave+s; diff=2*s;
	  }
	  crst[m]=(d.cn[m]-1)/diff;
	}
      }

      bool flag=true;
      for(map<unsigned short, short>::iterator i=num.begin(); i!=num.end(); ++i){
	line & s = d.sc[i->second];
	int n[2];
	for(int m=0; m<2; m++){
	  n[m]=lroundf((ctr(s, m)-cl[m])*crst[m]);
	  if(n[m]<0 && n[m]>=d.cn[m]){ cerr<<"Error in cell initialization"<<endl; exit(1); }

	  float d1=fabs(ctr(s, m)-(cl[m]+(n[m]-0.5f)/crst[m]));
	  float d2=fabs(ctr(s, m)-(cl[m]+(n[m]+0.5f)/crst[m]));
	  float d=min(d1, d2)*sin12-s.r;
	  if(d<0){ flag=false; cerr<<"Warning: string "<<(int)i->first<<" too close to cell boundary"<<endl; }
	}

	cells[n[0]][n[1]][i->first]++;
      }
      if(flag) d.rx=0;

      for(int m=0; m<2; m++){ d.cl[m]=cl[m]; d.crst[m]=crst[m]; }
    }

    {
      unsigned int pos=0;
      for(int i=0; i<d.cn[0]; i++) for(int j=0; j<d.cn[1]; j++){
	map<unsigned short, int> & c = cells[i][j];

	if(c.size()>0){
	  d.is[i][j]=pos;
	  for(map<unsigned short, int>::const_iterator n=c.begin(); n!=c.end(); ++n){
	    if(pos==NSTR){ cerr<<"Number of string cells exceeds capacity of "<<NSTR<<endl; exit(1); }
	    d.ls[pos++]=num[n->first];
	  }
	  d.ls[pos-1]|=0x8000;
	}
	else d.is[i][j]=0x8000;
      }
    }

    cerr<<"Loaded "<<d.gsize<<" DOMs ("<<d.cn[0]<<"x"<<d.cn[1]<<")"<<endl;

    d.tmod=0;

    { // initialize 2d ice layer tilt
      ifstream inFile((tiltdir+"tilt.set").c_str(), ifstream::in);
      if(!inFile.fail()){
	float mdir[2];
	int n=0;
	string in;
	while(getline(inFile, in)){
	  int tot=sscanf(in.c_str(), "%f %f %f %d %d", &mdir[n], &d.mmin[n], &d.mstp[n], &d.mnum[n], &d.mcut[n]);
	  if(tot==5) n++; if(n>=2) break;
	}
	inFile.close();
	if(n<2){ cerr << "File tilt.set found, but is corrupt" << endl; exit(1); }

	if(mdir[0]!=DIR1 || mdir[1]!=DIR2){ cerr << "Unsupported tilt grid configuration" << endl; exit(1); }

	int size=0;
	for(int j=0; j<d.mnum[1]; j++) for(int i=0; i<d.mnum[0]; i++) if(i<j+d.mcut[0] && j<i+d.mcut[1]) size++;
	if(size>LMAX){ cerr << "File tilt.set defines too many dust maps" << endl; exit(1); }
	if(d.mnum[0]>CTX || d.mnum[1]>CTY){ cerr<<"Error: tilt map conifguration exceeds buffer"<<endl; exit(1); }
	if(d.mstp[0]<=0 || d.mstp[1]<=0){ cerr << "Tilt map does not use increasing range order" << endl; exit(1); }

	for(int m=0, j=0; j<d.mnum[1]; j++) for(int i=0; i<d.mnum[0]; i++){
	    if(i<j+d.mcut[0] && j<i+d.mcut[1]) d.mcol[i][j]=m++;
	    else d.mcol[i][j]=-1;
	  }

	ifstream inFile((tiltdir+"tilt.map").c_str(), ifstream::in);
	if(!inFile.fail()){
	  d.lnum=size;
	  vector<float> pts(d.lnum), ds;
	  vector<float> lp[d.lnum];

	  while(getline(inFile, in)){
	    stringstream str(in);
	    float depth;
	    while(str >> depth){
	      int i=0;
	      while(str >> pts[i++]) if(i>=d.lnum) break;
	      if(i!=d.lnum) break;
	      ds.push_back(depth);
	      for(i=0; i<d.lnum; i++) lp[i].push_back(pts[i]);
	    }
	  }
	  inFile.close();

	  int size=ds.size();
	  if(size-1>LYRS){ cerr << "File tilt.map defines too many map points" << endl; exit(1); }
	  for(int i=1; i<size; i++) if(ds[i]<ds[i-1]){ cerr << "Tilt map does not use increasing depth order" << endl; exit(1); }
	  for(int i=0; i<d.lnum; i++) for(int j=0; j<size-1; j++){
	      z.lp[i][j].x=lp[i][size-1-j]; z.lp[i][j].y=lp[i][size-2-j]-lp[i][size-1-j];
	    }
	  d.lpts=size-2;

	  if(size<2) d.lnum=0;
	  else{
	    float lmin=ds[0], lmax=ds[size-1];
	    d.lmin=zoff-lmax; d.lrdz=(size-1)/(lmax-lmin);
	  }

	  if(d.lnum>0){
	    d.tmod=2;
	    cerr<<"Loaded "<<d.lnum<<"x"<<size<<" 2d dust layer points"<<endl;
	  }
	}
      }
    }

    if(d.lnum<=0){ // initialize 1d ice layer tilt
      d.lnum=0; d.l0=0, d.r0=0;
      const float thx=225;
      d.lnx=cos(fcv*thx), d.lny=sin(fcv*thx);

      ifstream inFile((tiltdir+"tilt.par").c_str(), ifstream::in);
      if(!inFile.fail()){
	int str;
	float aux;
	vector<float> lr;
	while(inFile >> str >> aux){ if(aux==0) d.l0=str; lr.push_back(aux); }
	inFile.close();

	int size=lr.size();
	if(size>LMAX){ cerr << "File tilt.par defines too many dust maps" << endl; exit(1); }
	for(int i=1; i<size; i++) if(lr[i]<lr[i-1]){ cerr << "Tilt map does not use increasing range order" << endl; exit(1); }
	for(int i=0; i<size; i++) d.lr[i]=lr[i];

	string in;
	ifstream inFile((tiltdir+"tilt.dat").c_str(), ifstream::in);
	if(!inFile.fail()){
	  d.lnum=size;
	  vector<float> pts(d.lnum), ds;
	  vector<float> lp[d.lnum];

	  while(getline(inFile, in)){
	    stringstream str(in);
	    float depth;
	    while(str >> depth){
	      int i=0;
	      while(str >> pts[i++]) if(i>=d.lnum) break;
	      if(i!=d.lnum) break;
	      ds.push_back(depth);
	      for(i=0; i<d.lnum; i++) lp[i].push_back(pts[i]);
	    }
	  }
	  inFile.close();

	  int size=ds.size();
	  if(size-1>LYRS){ cerr << "File tilt.dat defines too many map points" << endl; exit(1); }
	  for(int i=1; i<size; i++) if(ds[i]<ds[i-1]){ cerr << "Tilt map does not use increasing depth order" << endl; exit(1); }
	  for(int i=0; i<d.lnum; i++) for(int j=0; j<size-1; j++){
	      z.lp[i][j].x=lp[i][size-1-j]; z.lp[i][j].y=lp[i][size-2-j]-lp[i][size-1-j];
	    }
	  d.lpts=size-2;

	  if(size<2) d.lnum=0;
	  else{
	    float lmin=ds[0], lmax=ds[size-1];
	    d.lmin=zoff-lmax; d.lrdz=(size-1)/(lmax-lmin);
	  }

	  if(d.lnum>0){
	    d.tmod=1;
	    cerr<<"Loaded "<<d.lnum<<"x"<<size<<" dust layer points"<<endl;
	  }
	}
      }
    }

    { // should follow both tilt and geometry initializations
      int nr=0;
      float r0=0;
      for(int n=0; n<d.gsize; n++){
	int str=q.names[n].str;
	if(str==d.l0){
	  nr++;
	  DOM & om = q.oms[n];
	  r0+=d.lnx*om.r[0]+d.lny*om.r[1];
	}
      }
      if(nr>0) d.r0=r0/nr;
    }

    {
      char * VTHK=getenv("VTHK");
      d.vthk=VTHK==NULL?0:atoi(VTHK);
      cerr<<"Ice layer thickness: ";
      switch(d.vthk){
      case 0:
	cerr<<"UNIFORM";
	break;
      case 1:
	cerr<<"TILT-CORRECTED";
	break;
      default:
	cerr<<"UNKNOWN";
      }
      cerr<<endl;
    }

    { // initialize ice parameters
      int size;            // size of kurt table
      float dh, rdh, hmin; // step, 1/step, and min depth

      {
	vector<float> wx, wy, wz;
	bool flag=true, rdef=false;

	char * WFLA=getenv("WFLA");
	float wfla=WFLA==NULL?0:atof(WFLA);
	if(wfla>0){
	  for(float x=0.f; x<1.f+XXX; x+=0.1f){
	    wx.push_back(3*x*x-2*x*x*x);
	    wy.push_back(wfla+10.f*(2*x-1)); // +- 10 nm smoothstep
	  }
	  cerr<<"Using single wavelength="<<wfla<<" [nm]"<<endl;
	  qwv.ssw(wfla);
	}
	else{
	  string flwl("dat");
	  {
	    char * env=getenv("FLWL");
	    if(env!=NULL){
	      flwl=string(env);
	      cerr<<"wv.dat file extension is modified with FLWL="<<flwl<<endl;
	    }
	  }

	  flwl="wv."+flwl;
	  ifstream inFile((ppcdir+flwl).c_str(), ifstream::in);
	  if(!inFile.fail()){
	    int num=0;
	    float xa, ya, xo=0, yo=0;
	    while(inFile>>xa>>ya){
	      if(( xa<0 || 1<xa ) || (num==0 && xa!=0) || (num>0 && ( xa<=xo || ya<=yo ))){ flag=false; break; }
	      wx.push_back(xa); wy.push_back(ya);
	      xo=xa; yo=ya; num++;
	    }
	    if(xo!=1 || wx.size()<2) flag=false;
	    inFile.close();
	    if(flag){ cerr<<"Loaded "<<wx.size()<<" wavelenth points"<<endl; }
	    else{ cerr<<"File "<<flwl<<" did not contain valid data"<<endl; exit(1); }
	  }
	  else{ cerr<<"Could not open file "<<flwl<<endl; exit(1); }
	}

	if(flag){
	  vector<float> qx, qy;

	  ifstream inFile((ppcdir+"wv.rde").c_str(), ifstream::in);
	  if(!inFile.fail()){
	    int num=0;
	    bool flag=true;
	    float xa, ya, yo=0;
	    while(inFile>>ya>>xa){
	      if(xa<0 || (num>0 && ya<=yo)){ flag=false; break; }
	      qx.push_back(xa); qy.push_back(ya);
	      yo=ya; num++;
	    }
	    if(qx.size()<2) flag=false;
	    inFile.close();
	    if(flag){ cerr<<"Loaded "<<qx.size()<<" RDE coefficients"<<endl; }
	    else{ cerr<<"File wv.rde did not contain valid data"<<endl; exit(1); }
	    rdef=flag;
	  }

	  if(rdef){
	    int k=0, n=qy.size();
	    for(vector<float>::iterator i=wy.begin(); i!=wy.end(); i++){
	      float w=*i, r;
	      for(; k<n; k++) if(qy[k]>w) break;
	      if(k==0) r=qx[0];
	      else if(k==n) r=qx[n-1];
	      else{
		r=((w-qy[k-1])*qx[k]+(qy[k]-w)*qx[k-1])/(qy[k]-qy[k-1]);
		k--;
	      }
	      wz.push_back(r);
	    }
	  }

	  {
	    vector<float> qw, qf, qz;
	    float qo=0, wo=wy[0];
	    unsigned int n=wx.size();
	    for(unsigned int i=1; i<n; i++){
	      float qv=wx[i], wv=wy[i];
	      float bin=wv-wo;
	      float wva=wv-bin/2;
	      float val=dppm*doma*omav*(qv-qo)/bin;
	      qw.push_back(wva);
	      qf.push_back(val);
	      if(rdef) qz.push_back(val*(wz[i]+wz[i-1])/2);
	      qo=qv, wo=wv;
	    }
	    float blo=wy[1]-wy[0], bhi=wy[n-1]-wy[n-2];
	    irdes[make_pair(-1,0)].read(qw, qf, blo, bhi);
	    if(rdef) irdes[make_pair(-1,1)].read(qw, qz, blo, bhi);
	  }
	}
      }

      for(map<pair<int,int>, irde>::iterator i=irdes.begin(); i!=irdes.end(); ++i){
	const pair<int,int> & m = i->first;
	if(m.first>=0){
	  stringstream str;
	  str<<omdir<<"om.wv_"<<m.first<<"."<<m.second;
	  const string & file=str.str();

	  bool flag=false;
	  vector<float> qw, qf;

	  ifstream inFile(file.c_str(), ifstream::in);
	  if(!inFile.fail()){
	    int num=0;
	    flag=true;
	    float w, f, u=0;
	    while(inFile>>w>>f){
	      if(f<0 || (num>0 && w<=u)){ flag=false; break; }
	      qf.push_back(f*qwv.cherenkov(w)); qw.push_back(w);
	      u=w; num++;
	    }
	    if(qf.size()<2) flag=false;
	    inFile.close();
	    if(flag){ cerr<<"Loaded "<<qf.size()<<" wv points from file "<<file<<endl; }
	    else{ cerr<<"File "<<file<<" did not contain valid data"<<endl; exit(1); }
	  }

	  if(flag) i->second.read(qw, qf);
	  else cerr<<"Warning: file "<<file<<" was not read. Omitting hits in this OM type!"<<endl;
	}
      }

      {
	float sum=0;
	int k=0;
	for(float w=qwv.wmin+qwv.bin/2; w<qwv.wmax; w+=qwv.bin, k++){
	  float val=0;
	  for(map<pair<int,int>, irde>::iterator i=irdes.begin(); i!=irdes.end(); ++i){
	    irde & mt = i->second;
	    if(!mt.dat.empty()){
	      float ri=mt.dat[k]*types[i->first.first].rde*mt.rmax*qwv.bin;
	      mt.rde+=ri; mt.sum.push_back(mt.rde); mt.dat[k]=ri;
	      if(val<ri) val=ri;
	    }
	  }
	  sum+=val;
	  env.dat.push_back(val);
	  env.sum.push_back(sum);
	}
	env.rde=sum; q.eff*=sum;
	cerr<<"Calculated value of ppm="<<q.eff<<endl;

	for(int i=0; i<k; i++) env.sum[i]/=sum;
      }

      float wv0=400;
      float A, B, D, E, a, k;
      float Ae, Be, De, Ee, ae, ke;
      vector<float> dp, be, ba, td, k1, k2, ra, rb;

      {
	bool flag=true, fail=false;
	ifstream inFile((icedir+"icemodel.par").c_str(), ifstream::in);
	if((flag=!inFile.fail())){
	  if(flag) flag=(bool)(inFile >> a >> ae);
	  if(flag) flag=(bool)(inFile >> k >> ke);
	  if(flag) flag=(bool)(inFile >> A >> Ae);
	  if(flag) flag=(bool)(inFile >> B >> Be); fail=!flag;
	  if(flag) flag=(bool)(inFile >> D >> De); D=flag?D/pow(wv0, k):1;
	  if(flag) flag=(bool)(inFile >> E >> Ee); E=flag?E/pow(wv0, k):0;
	  if(fail) cerr << "File icemodel.par found, but is corrupt" << endl;
	  inFile.close(); if(fail) exit(1);
	}
	else{ cerr << "File icemodel.par was not found" << endl; exit(1); }
      }

      {
	ifstream inFile((icedir+"icemodel.dat").c_str(), ifstream::in);
	if(!inFile.fail()){
	  size=0;
	  float dpa, bea, baa, tda, k1a, k2a, bfra, bfrb;

	  string in;
	  while(getline(inFile, in)){
	    int num=sscanf(in.c_str(), "%f %f %f %f %f %f %f %f", &dpa, &bea, &baa, &tda, &k1a, &k2a, &bfra, &bfrb);

	    if(num>=4){
	      dp.push_back(dpa);
	      be.push_back(bea);
	      ba.push_back(baa);
	      td.push_back(tda);
	      k1.push_back(num>=5?exp(k1a):dk1);
	      k2.push_back(num>=6?exp(k2a):dk2);
	      ra.push_back(num>=7?bfra:1);
	      rb.push_back(num>=8?bfrb:1);
	      size++;
	    }
	  }
	  inFile.close();
	  if(size<1){ cerr << "File icemodel.dat found, but is corrupt" << endl; exit(1); }
	}
	else{ cerr << "File icemodel.dat was not found" << endl; exit(1); }
      }

      dh=size>1?(dp[size-1]-dp[0])/(size-1):100.f;
      if(dh<=0){ cerr << "Ice table does not use increasing depth spacing" << endl; exit(1); }

      for(int i=0; i<size; i++) if(i>0) if(fabs(dp[i]-dp[i-1]-dh)>dh*XXX){
	cerr << "Ice table does not use uniform depth spacing" << endl; exit(1);
      }
      cerr<<"Loaded "<<size<<" ice layers"<<endl;

      if(size>MAXLYS){
	cerr<<"Error: too many layers ("<<size<<"), truncating to "<<MAXLYS<<endl;
	size=MAXLYS;
      }

      rdh=1/dh; hmin=zoff-dp[size-1];
      {
	d.size=size;
	d.dh=dh;
	d.rdh=rdh;
	d.hmin=hmin;
      }

      float bble=0, bblz, bbly;

      {
	ifstream inFile((icedir+"icemodel.bbl").c_str(), ifstream::in);
	if(!inFile.fail()){
	  if(!(inFile >> bble >> bblz >> bbly)){
	    cerr << "File icemodel.bbl found, but is corrupt" << endl;
	    bble=0;
	  }
	  else{
	    cerr << "Air bubble parameters: " << bble << " " << bblz << " " << bbly << endl;
	  }
	}
      }

      for(int i=0; i<size; i++){
	int j=size-1-i;
	d.az[i].k1=k1[j]; d.az[i].k2=k2[j];
	d.az[i].ra=ra[j]; d.az[i].rb=ra[j]*rb[j];
      }

      float arf=0;
      float srw=0, srf=1;

      {
	int version=0;
	ifstream inFile((icedir+"version").c_str(), ifstream::in);
	if(!inFile.fail()){
	  if(!(inFile >> version)){
	    cerr << "File version found, but is corrupt" << endl;
	    version=0;
	  }
	  else{
	    cerr << "Adjusting defaults for ice model version = " << version << endl;
	  }
	}

	switch(version){
	case 10: arf=0; srw=0, srf=0; break;
	case 20: arf=1; srw=1, srf=0; break;
	}
      }

      {
	char * absm=getenv("ABSM");
	if(absm!=NULL){
	  cerr<<"Ice absorption table: ";
	  switch(atoi(absm)){
	  case 0:
	    arf=0;
	    cerr<<"MieDUST";
	    break;
	  case 1:
	    arf=1;
	    cerr<<"FULL400";
	    break;
	  default:
	    cerr<<"UNKNOWN";
	  }
	  cerr<<endl;
	}
      }

      {
	char * bfrm=getenv("BFRM");
	if(bfrm!=NULL){
	  cerr<<"BFR scattering subtraction method: ";
	  switch(atoi(bfrm)){
	  case 0:
	    srw=0, srf=1;
	    cerr<<"DEFAULT";
	    break;
	  case 1:
	    srw=1, srf=0;
	    cerr<<"CORRECT";
	    break;
	  case 2:
	    srw=0, srf=0;
	    cerr<<"NO SUB.";
	    break;
	  default:
	    cerr<<"UNKNOWN";
	  }
	  cerr<<endl;
	}
      }

      for(int n=0; n<WNUM; n++){
	float wvi=env.binf((n+0.f)/WNUM);
	int mi=min(max((int) floor(wvi), 0), qwv.num-1);

	float wvf=env.binf((n+1.f)/WNUM);
	int mf=min(max((int) floor(wvf), 0), qwv.num-1);

	float wva=qwv.wav(env.binf((n+0.5f)/WNUM));
	q.wvs[n].w=wva; q.wvs[n].i=qwv.wav(wvi), q.wvs[n].f=qwv.wav(wvf);

	for(map<pair<int,int>, irde>::iterator i=irdes.begin(); i!=irdes.end(); ++i){
	  irde & mt = i->second;
	  if(mt.rde>0){
	    vector<float> & wy = mt.dat;
	    float rde;
	    if(mf==mi) rde=wy[mi]*(wvf-wvi);
	    else{
	      rde=(mi+1-wvi)*wy[mi]+(wvf-mf)*wy[mf];
	      for(int i=mi+1; i<mf; i++) rde+=wy[i];
	    }

	    rde/=env.rde/WNUM;
	    mt.rat.push_back(rde);
	  }
	}

	float l_a=pow(wva/wv0, -a);
	float l_k=pow(wva/wv0, -k);
	float AB0=A*exp(-B/wv0);
	float ABl=A*exp(-B/wva);

	ices & w = z.w[n];

	for(int i=0; i<size; i++){
	  int j=size-1-i;
	  float bbl=bble>0?dp[j]<bblz&&dp[j]<bbly?bble*(1-dp[j]/bblz)*(1-dp[j]/bbly):0:0;
	  float sca = (bbl+be[j]*l_a - ra[j]*d.sum*(srf+srw*l_a))/(1-d.g);
	  float abs = (D*ba[j]+E)*l_k + (ABl-arf*AB0*l_k)*(1+0.01*td[j]);
	  if(sca>0 && abs>0) w.z[i].sca=sca, w.z[i].abs=abs;
	  else{ cerr << "Invalid value of ice parameter, cannot proceed" << endl; exit(1); }
	}

	float ng, np=qwv.np(wva*1.e-3, &ng);
	float c=0.299792458; d.ocv=1/c; w.wvl=n; w.ocm=ng/c;
	w.coschr=1/np; w.sinchr=sqrt(1-w.coschr*w.coschr);
      }
    }

    cerr<<endl;
  }
} m;
