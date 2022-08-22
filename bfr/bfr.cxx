#include <cmath>
#include <cstdlib>
#include <vector>

#include <fstream>
#include <iostream>

#include <gsl/gsl_multifit.h>

using namespace std;

bool verbose=false;
bool girdle=true;
bool orefr=false;
bool orefl=false;
bool loop=true; // guarantee interaction at every step
double frac=0.0;
double elong=1.0;
double xx=1.e-10;

double xrnd(){
  return drand48();
}

class myvect{
public:
  double x, y, z;

  myvect() : x(0), y(0), z(0) {}

  myvect(double x, double y, double z){
    this->x=x, this->y=y, this->z=z;
  }

  void set(myvect q){
    x=q.x, y=q.y, z=q.z;
  }

  void divide(double p){
    x/=p, y/=p, z/=p;
  }

  void multiply(double p){
    x*=p, y*=p, z*=p;
  }

  myvect operator* (double p){
    return myvect(x*p, y*p, z*p);
  }

  myvect operator/ (double p){
    return myvect(x/p, y/p, z/p);
  }

  myvect operator+ (myvect q){
    return myvect(x+q.x, y+q.y, z+q.z);
  }

  myvect operator- (myvect q){
    return myvect(x-q.x, y-q.y, z-q.z);
  }

  myvect smallest(){
    myvect temp(fabs(x), fabs(y), fabs(z));
    switch( temp.z<temp.y? temp.z<temp.x?2:0 : temp.y<temp.x?1:0 ){
    case 0: temp=myvect(1, 0, 0); break;
    case 1: temp=myvect(0, 1, 0); break;
    case 2: temp=myvect(0, 0, 1); break;
    }
    return temp;
  }

  myvect cross(myvect q){
    return myvect(y*q.z-q.y*z, z*q.x-q.z*x, x*q.y-q.x*y);
  }

  double dot(myvect q){
    return x*q.x+y*q.y+z*q.z;
  }

  double norm(){
    return sqrt(this->dot(*this));
  }

  void normalize(){
    divide(norm());
  }

  myvect np(myvect q){ // normalized vector perpendicular to this and q (assume |q|=1)
    // q.normalize();
    myvect temp=this->cross(q);
    double size=temp.norm();
    if(size>0) temp.divide(size);
    else{
      temp=this->smallest().cross(q);
      temp.normalize();
    }
    return temp;
  }
};

ostream& operator<<(ostream& o, const myvect& q){
  o << "(" << q.x << "," << q.y << "," << q.z << ")";
  return o;
}

class surface:
  public myvect // normal vector pointing first into second medium
{
public:
  surface() : myvect() {}
  surface(myvect q) : myvect(q) {}
  surface(double x, double y, double z) : myvect(x, y, z) {}

  myvect p1, p2, e;
  bool skip; // no crossing this iteration (skip this step)

  void setp(myvect q){ // (assume |*this|=1)
    q.normalize();
    p1=this->np(q);
    p2=this->cross(p1);
    // p2.normalize();
  }

  void setp(){ // (assume |*this|=1)
    p1=this->cross(this->smallest());
    p1.normalize();
    p2=this->cross(p1);
    // p2.normalize();
  }

  // random surface segment of an ellipsoid, representing the average grain
  // e.x,y,z=a,b,c are the x,y,z radia
  // computes gradients of ellipsoid surface from x^2/a^2+y^2/b^2+z^2/c^2=1 defining equation
  // for directions sampled from sphere and weights according to
  // https://math.stackexchange.com/questions/973101/how-to-generate-points-uniformly-distributed-on-the-surface-of-an-ellipsoid
  myvect ellipsoid(){  // substitute for rand_x().
    double weight = 0;
    double maxweight = min(e.x, min(e.y, e.z));
    myvect r;

    do {
      double costheta=2*xrnd()-1;
      double sintheta=sqrt(1-costheta*costheta);
      double ph=2*M_PI*xrnd();
      double cosphi=cos(ph), sinphi=sin(ph);

      r.x=sintheta*cosphi/e.x;
      r.y=sintheta*sinphi/e.y;
      r.z=costheta/e.z;

      weight=maxweight*r.norm();
    } while((skip=(xrnd() >= weight)) && loop);

    r.normalize();
    r=p1*r.x+p2*r.y+(*this)*r.z;
    return r;
  }

  myvect rand(){ // random from sphere
    double ct=2*xrnd()-1;
    double st=sqrt(1-ct*ct);
    double ph=2*M_PI*xrnd();
    double cp=cos(ph), sp=sin(ph);
    return myvect(st*cp, st*sp, ct);
  }

  double elong_sampling(double p){ // returns the cos(th) value for elongated (stretched) ice
    double p2=p*p, pm=max(1., p);
    double weight, xi, area;

    do {
      xi=2*xrnd()-1;
      area=sqrt(p2+(1-p2)*xi*xi);
      weight=area/pm;
    } while((skip=(xrnd() >= weight)) && loop);

    return xi/area;
  }

  myvect rand_x(){ // interface plane orientation
    if(elong<0) return ellipsoid();
    myvect r;
    if(elong!=1){
      double ct=elong_sampling(elong);
      double st=sqrt(1-ct*ct);
      double ph=2*M_PI*xrnd();
      double cp=cos(ph), sp=sin(ph);
      r=(*this)*ct+(p1*cp+p2*sp)*st;
    }
    else r=rand(), skip=false;
    return r;
  }

  myvect rand_i(myvect q){
    // gets grain boundary plane according to uniform distribution on sphere (ellipsoid optional)
    // dot product between boundary and poynting vector = probability to see plane
    myvect r;
    do r=rand_x(); while((!skip) && (skip=(xrnd() >= fabs(r.dot(q)/q.norm()))) && loop);
    return r;
  }

  myvect r_sav;

  myvect rand_c(bool reset=false){ // c-axis oriention
    myvect r;
    if(frac>0){
      if(reset){
	double ph=2*M_PI*xrnd();
	double cp=cos(ph), sp=sin(ph);
	r_sav=p1*cp+p2*sp;
      }
      r=xrnd()<frac?r_sav:rand();
    }
    else{
      if(girdle){
	double ph=2*M_PI*xrnd();
	double cp=cos(ph), sp=sin(ph);
	r=p1*cp+p2*sp;
      }
      else r=rand();
    }
    return r;
  }
};

class photon:
  public myvect // coordinates
{
public:
  photon() : myvect() {}
  photon(myvect q) : myvect(q) {}
  photon(double x, double y, double z) : myvect(x, y, z) {}

  double n;
  myvect s; // photon propagation direction
  myvect k, ki; // wave vector
  myvect H, Hi; // magnetic field
  myvect E, Ei; // electric field
  myvect D, Di; // displacement field and polarization

  void advance(double r){
    set(* this + s * (r/s.norm()));
  }
};

ostream& operator<<(ostream& o, const photon& q){
  o << "n = " << q.n << endl;
  o << "r = " << (myvect) q << endl;
  o << "s = " << q.s << endl;
  o << "k = " << q.k << " " << "ki = " << q.ki << endl;
  o << "H = " << q.H << " " << "Hi = " << q.Hi << endl;
  o << "E = " << q.E << " " << "Ei = " << q.Ei << endl;
  o << "D = " << q.D << " " << "Di = " << q.Di << endl;
  return o;
}

surface plane_saved;
class medium:
  public myvect // optical axis
{
public:
  double no,  ne; // refractive indices
  double beta;

  enum type{
    any = 3,
    ordinary = 1,
    extraordinary = 2
  } current;

  medium() : myvect() {
    set_n();
  }

  medium(myvect q) : myvect(q) {
    set_n();
  }

  medium(double x, double y, double z) : myvect(x, y, z) {
    set_n();
  }

  void set_n(){
    no=1.3185, ne=1.3200; // at 405 nm
    double aux=ne/no;
    beta=aux*aux-1;
    current=any;
  }

  void set_n(double n){
    no=n, ne=n, beta=0;
  }

  void eires(myvect q, myvect E, myvect D){
    double k2=q.dot(q);
    double nE=q.dot(E);
    myvect a=E*k2-q*nE-D;
    cout<<"EIGEN RESIDUAL NORM: "<<a.norm()<<endl;
  }

  void eires(myvect q, myvect qi, myvect E, myvect Ei, myvect D, myvect Di){
    double k2r=q.dot(q)-qi.dot(qi), k2i=2*q.dot(qi);
    double nEr=q.dot(E)-qi.dot(Ei);
    double nEi=qi.dot(E)+q.dot(Ei);
    myvect ar=(E*k2r-Ei*k2i)-(q*nEr-qi*nEi)-D;
    myvect ai=(E*k2i+Ei*k2r)-(qi*nEr+q*nEi)-Di;
    cout<<"EIGEN RESIDUAL NORM: "<<ar.norm()<<" "<<ai.norm()<<endl;
  }

  void set_k(myvect q, photon & o, photon & e, bool p = false){ // o: ordinary; e: extraordinary; q=p?s:k
    myvect& axis = *this;

    q.normalize();

    o.D=axis.np(q);
    o.n=no;
    o.k=q;
    o.k.multiply(o.n);

    if(p){
      myvect E=q.cross(o.D);
      myvect axe=axis*axis.dot(E);
      e.D=axe*(ne*ne)+(E-axe)*(no*no);
      e.D.normalize();
      q=o.D.cross(e.D);
    }
    else e.D=q.cross(o.D);

    double cs=axis.dot(q)/no, sn=axis.dot(e.D)/ne;
    e.n=1/sqrt(cs*cs+sn*sn);
    e.k=q;
    e.k.multiply(e.n);

    o.E=o.D/(no*no);
    o.H=o.k.cross(o.E);
    o.s=o.E.cross(o.H);

    myvect axd=axis*axis.dot(e.D);
    e.E=axd/(ne*ne)+(e.D-axd)/(no*no);
    e.H=e.k.cross(e.E);
    e.s=e.E.cross(e.H);

    if(verbose){
      if(current&ordinary){
	eires(o.k, o.E, o.D);
	cout<<"\tk="<<o.k<<"\ts="<<o.s<<"\tdot="<<o.k.dot(o.s)*o.n*o.n<<endl;
      }
      if(current&extraordinary){
	eires(e.k, e.E, e.D);
	cout<<"\tk="<<e.k<<"\ts="<<e.s<<"\tdot="<<e.k.dot(e.s)*e.n*e.n<<endl;
      }
    }
  }

  void set_k(myvect q, myvect qi, photon & o, photon & e){ // o: ordinary; e: extraordinary
    myvect& axis = *this;

    myvect tr=axis.cross(q), ti=axis.cross(qi);
    double sr=tr.norm(), si=ti.norm();
    double size=sqrt(sr*sr+si*si);
    if(size==0){
      myvect temp=axis.smallest();
      tr=temp.cross(q), ti=temp.cross(qi);
      sr=tr.norm(), si=ti.norm();
      size=sqrt(sr*sr+si*si);
    }
    tr.divide(size), ti.divide(size);

    o.D=tr, o.Di=ti;
    o.n=no;
    o.k=q, o.ki=qi;

    e.D=q.cross(o.D)-qi.cross(o.Di), e.Di=qi.cross(o.D)+q.cross(o.Di);
    sr=e.D.norm(), si=e.Di.norm();
    size=sqrt(sr*sr+si*si);
    if(size>0) e.D.divide(size), e.Di.divide(size);

    sr=q.norm(), si=qi.norm();
    e.n=sqrt(sr*sr+si*si);
    e.k=q, e.ki=qi;

    o.E=o.D/(no*no), o.Ei=o.Di/(no*no);
    o.H=o.k.cross(o.E)-o.ki.cross(o.Ei), o.Hi=o.ki.cross(o.E)+o.k.cross(o.Ei);
    o.s=o.E.cross(o.H)+o.Ei.cross(o.Hi);

    myvect axd=axis*axis.dot(e.D), axdi=axis*axis.dot(e.Di);
    e.E=axd/(ne*ne)+(e.D-axd)/(no*no), e.Ei=axdi/(ne*ne)+(e.Di-axdi)/(no*no);
    e.H=e.k.cross(e.E)-e.ki.cross(e.Ei), e.Hi=e.ki.cross(e.E)+e.k.cross(e.Ei);
    e.s=e.E.cross(e.H)+e.Ei.cross(e.Hi);

    if(verbose){
      if(current&ordinary){
	eires(o.k, o.ki, o.E, o.Ei, o.D, o.Di);
	cout<<"\tk="<<o.k<<"\tki="<<o.ki<<"\tdot="<<plane_saved.dot(o.s)<<endl;
      }
      if(current&extraordinary){
	eires(e.k, e.ki, e.E, e.Ei, e.D, e.Di);
	cout<<"\tk="<<e.k<<"\tki="<<e.ki<<"\tdot="<<plane_saved.dot(e.s)<<endl;
      }
    }
  }

  vector<photon> set_k(myvect k, surface s, bool same){
    vector<photon> result;

    double ky=s.dot(k);
    myvect Y=s;
    if(ky>0){
      ky*=-1;
      Y.multiply(-1);
    }

    myvect X=k-Y*ky;
    double kx=X.norm();
    if(kx>0) X.divide(kx);

    bool evan=true;
    {
      current=ordinary;
      double D=no*no-kx*kx;
      if(D>=0){
	D=sqrt(D);
	double ry=same?D:-D;
	myvect K=X*kx+Y*ry;
	photon O, E;
	set_k(K, O, E);
	result.push_back(O);
      }
      else if(evan){
	D=sqrt(-D);
	double ry=same?D:-D;
	myvect K=X*kx, Ki=Y*ry;
	photon O, E;
	set_k(K, Ki, O, E);
	result.push_back(O);
      }
    }

    {
      current=extraordinary;
      double ax=this->dot(X);
      double ay=this->dot(Y);
      double D=ne*ne*(1+beta*ay*ay)-kx*kx*(1+beta*(ax*ax+ay*ay));
      if(D>=0){
	D=sqrt(D);
	double b=-beta*ax*ay*kx, a=1+beta*ay*ay;
	double r1=(b+D)/a, r2=(b-D)/a;
	double ry=same?r1:r2;
	myvect K=X*kx+Y*ry;
	photon O, E;
	set_k(K, O, E);
	result.push_back(E);
      }
      else if(evan){
	D=sqrt(-D);
	double b=-beta*ax*ay*kx, a=1+beta*ay*ay;

	myvect K=X*kx+Y*(b/a);
	double r1=D/a, r2=-D/a;
	double ry=same?r1:r2;
	myvect Ki=Y*ry;

	photon O, E;
	set_k(K, Ki, O, E);
	result.push_back(E);
      }
    }

    current=any;
    return result;
  }
};

bool interact(medium one, medium two, surface plane, photon & p){
  myvect r(p);
  // plane.setp(); //  nominal  basis
  plane.setp(p.k); // optimized basis

  if(verbose){
    cout<<"one: "<<one<<endl;
    cout<<"two: "<<two<<endl;
    cout<<"ref: "<<plane<<" "<<plane.p1<<" "<<plane.p2<<endl;
  }

  bool same=true;
  if(fabs(p.s.dot(plane))<xx){
    cerr<<"photon collinear to interface surface"<<endl;
    return same;
  }

  plane_saved=plane;
  vector<photon> reflected=one.set_k(p.k, plane, true);
  vector<photon> refracted=two.set_k(p.k, plane, false);

  int refl=reflected.size();
  int refr=refracted.size();
  int dim=refl+refr, num=6;
  int xdm=2*dim, xnm=2*num;
  double chisq=0;

  if(dim<=0){
    cerr<<"no dice:"<<endl<<"1 = "<<one<<endl<<"2 = "<<two<<endl<<"i = "<<plane<<endl<<p<<endl;
    return same;
  }

  double q[xnm][dim], f[xnm], x[dim], s[dim+1];

  for(int i=0; i<=dim; i++){
    for(int j=0; j<xnm; j++) f[j]=0;
    photon& tmp = i<refl ? reflected[i] : i<dim ? refracted[i-refl] : p;

    if(verbose) cout << i << ": " << tmp << endl;

    f[0]=tmp.D.dot(plane);
    f[1]=tmp.E.dot(plane.p1);
    f[2]=tmp.E.dot(plane.p2);

    f[3]=tmp.H.dot(plane);
    f[4]=tmp.H.dot(plane.p1);
    f[5]=tmp.H.dot(plane.p2);

    f[6]=tmp.Di.dot(plane);
    f[7]=tmp.Ei.dot(plane.p1);
    f[8]=tmp.Ei.dot(plane.p2);

    f[9]=tmp.Hi.dot(plane);
    f[10]=tmp.Hi.dot(plane.p1);
    f[11]=tmp.Hi.dot(plane.p2);

    if(i<dim) for(int j=0; j<xnm; j++) q[j][i]=(i<refl?-1:1)*f[j];
    s[i]=tmp.s.dot(plane); // norm();
  }

  {
    gsl_matrix * A = gsl_matrix_alloc(xnm, xdm);
    gsl_vector * B = gsl_vector_alloc(xnm);

    for(int i=0; i<xnm; i++){
      for(int j=0; j<xdm; j++) gsl_matrix_set(A, i, j, j<dim?q[i][j]:i<num?-q[i+num][j-dim]:q[i-num][j-dim]);
      gsl_vector_set(B, i, f[i]);
    }

    if(verbose) for(int i=0; i<num; i++) cout<<"f["<<i<<"]: "<<f[i]<<" "<<f[i+num]<<endl;

    gsl_multifit_linear_workspace * W = gsl_multifit_linear_alloc(xnm, xdm);
    gsl_vector * X = gsl_vector_alloc(xdm);
    gsl_matrix * C = gsl_matrix_alloc(xdm, xdm);

    gsl_multifit_linear(A, B, X, C, &chisq, W);

    /*
    size_t rank;
    double tol=1.e-8;
    gsl_multifit_linear_svd(A, B, tol, & rank, X, C, &chisq, W); cout<<"rank="<<rank<<endl;
    */

    for(int i=0; i<dim; i++){
      double xre=gsl_vector_get(X, i), xim=gsl_vector_get(X, i+dim);
      x[i]=xre*xre+xim*xim; if(verbose) cout<<"x["<<i<<"]: "<<xre<<" "<<xim<<endl;
    }

    gsl_vector_free(X);
    gsl_matrix_free(C);
    gsl_multifit_linear_free(W);

    gsl_vector_free(B);
    gsl_matrix_free(A);
  }

  double sum=0;
  for(int i=0; i<dim; i++){ s[i]=fabs(x[i]*s[i]/s[dim]); sum+=s[i]; }
  if(verbose) cout<<refl<<" "<<refr<<"  "<<chisq<<" "<<sum<<" "<<p.k.dot(plane)<<" "<<p.s.dot(plane)<<endl;
  // cout<<chisq<<" "<<sum<<endl;

  if(verbose){
    cout<<"x:"; for(int i=0; i<dim; i++) cout<<" "<<x[i]; cout<<endl;
    cout<<"f:"; for(int i=0; i<dim; i++) cout<<" "<<s[i]; cout<<endl;
  }

  if(orefr){
    sum=0;
    for(int i=0; i<refl; i++) s[i]=0;
    for(int i=refl; i<dim; i++) sum+=s[i];
  }

  if(orefl){
    sum=0;
    for(int i=0; i<refl; i++) sum+=s[i];
    for(int i=refl; i<dim; i++) s[i]=0;
  }

  double y=sum*xrnd();
  sum=0;

  for(int i=0; i<dim; i++){
    sum+=s[i];
    if(y<sum){
      if(i<refl){
	same=true;
	p=reflected[i];
      }
      else{
	same=false;
	p=refracted[i-refl];
      }
      break;
    }
  }

  if(verbose) cout<<endl<<endl;
  p.set(r);
  return same;
}

void test(){
  verbose=true;

  myvect k(0, 0, 1);
  medium one(1, 0, 0), two(0, 1, 0); // two media definitions
  double f=0.03;
  surface plane(sqrt(1-f*f), 0, f);  // interface between two media

  /*
    myvect k(0.0953507,0.00228455,1.31646); k.normalize();
    medium one(0.260065,0.92991,-0.260065), two(-0.601485,0.525768,0.601485);
    surface plane(0.927202,0.371614,-0.0469062);
    one.normalize();
    two.normalize();
    plane.normalize();
  */

  // medium one(0, 0, 1), two(0, 0, 1);
  // one.set_n(2), two.set_n(1);

  photon o, e;
  one.set_k(k, o, e, true);
  if(verbose){
    cout << "o:" << endl << o << endl;
    cout << "e:" << endl << e << endl;
  }

  for(int n=0; n<2; n++){
    photon p=n==0?o:e; // incident photon
    interact(one, two, plane, p);
  }
}

void test2(double p, int num, int tot){
  // verbose=true;

  myvect ki(0, 0, 1);
  surface flow(sin(p), 0, cos(p)); // direction of ice flow
  flow.setp();
  surface crys=flow;

  if(elong<0){
    myvect pq[3];
    ifstream inFile("/dev/stdin", ifstream::in);
    if(!inFile.fail()){
      int size=0;
      double x, y, z;

      string in;
      while(getline(inFile, in)){
	int num=sscanf(in.c_str(), "%lf %lf %lf", &x, &y, &z);
	if(num==3) pq[size++]=myvect(x, y, z);
      }
      inFile.close();
      if(size!=3){
	cerr << "expecting the following on stdin:" << endl;
	cerr << "a b c [ ellipsoid ]" << endl;
	cerr << "x y z [ ice  flow ]" << endl;
	cerr << "x y z [photon sdir]" << endl;
	exit(1);
      }
    }
    else{ cerr << "cannot read from stdin" << endl; exit(1); }

    {
      crys=myvect(0, 0, 1);
      crys.p1=myvect(1, 0, 0);
      crys.p2=myvect(0, 1, 0);
      crys.e=pq[0];
      flow=pq[1];
      flow.setp();
      ki=pq[2];
    }
  }

  for(int j=0; j<tot; j++){
    myvect k(ki);
    medium one, two; // two media definitions
    surface plane; // interface between two media

    one = flow.rand_c(true);
    myvect r(0, 0, 0);
    photon o(r), e(r);
    one.set_k(k, o, e, true);
    photon p=xrnd()<0.5?o:e;

    if(verbose) cout<<endl<<endl;
    two = flow.rand_c();

    int count=0;
    for(int i=0; i<num; i++){
      p.advance(1./num);
      plane = crys.rand_i(p.s);
      if(crys.skip) continue;
      count++;
      bool same=interact(one, two, plane, p);
      if(!same) one=two;
      two = flow.rand_c();
    }

    myvect s=p.s;
    s.normalize();
    cout<<s.x<<" "<<s.y<<" "<<s.z<<" "<<p.x<<" "<<p.y<<" "<<p.z<<" "<<count<<endl;
  }
}

int main(int arg_c, char *arg_a[]){
  {
    char * bmp=getenv("SEED");
    if(bmp!=NULL){
      int seed=atoi(bmp);
      srand48(seed);
      cerr<<"SEED="<<seed<<endl;
    }
  }
  {
    char * bmp=getenv("ORFR");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)){
      orefr=true;
      cerr<<"Only refraction; no reflection"<<endl;
    }
  }
  {
    char * bmp=getenv("ORFL");
    if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)){
      orefl=true;
      cerr<<"Only reflection; no refraction"<<endl;
    }
  }

  if(arg_c>1){
    {
      char * bmp=getenv("NOLO");
      if(bmp!=NULL && !(*bmp!=0 && atoi(bmp)==0)){
	loop=false;
	cerr<<"Do not force crossings"<<endl;
      }
    }
    int inum=1000;
    {
      char * bmp=getenv("INUM");
      if(bmp!=NULL) inum=atoi(bmp);
      cerr<<"INUM="<<inum<<endl;
    }
    int jnum=100000;
    {
      char * bmp=getenv("JNUM");
      if(bmp!=NULL) jnum=atoi(bmp);
      cerr<<"JNUM="<<jnum<<endl;
    }

    if(arg_c>2) girdle = atof(arg_a[2])>0;
    if(arg_c>3)  elong = atof(arg_a[3]);
    if(arg_c>4)   frac = atof(arg_a[4]);

    test2(atof(arg_a[1])*M_PI/180, inum, jnum);
  }
  else{
    cerr<<"Usage: [SEED=0] [ORFR=0] [NOLO=0] [INUM=1000] [JNUM=100000] bfr [angle to flow] [girdle] [elongation] [unimodal]"<<endl;
    test();
  }
}
