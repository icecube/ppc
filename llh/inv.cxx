#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <map>
#include <vector>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

#include "../llh/nnls.cxx"
#include "../llh/lssl.cxx"

const int dim=1; //40
const double noise=25*500e-9;

typedef pair<int, int> key;

struct xyz{
  double x, y, z;

  xyz(){ x=0, y=0, z=0; }

  xyz(double x, double y, double z){
    this->x=x, this->y=y, this->z=z;
  }
};

main(int arg_c, char *arg_a[]){
  int srep=10, drep=240;
  double escl=1.0;

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
    char * bmp=getenv("ESCL");
    if(bmp!=NULL) escl=atof(bmp);
    cerr<<"ESCL="<<escl<<endl;
  }

  double dsr=drep/(double)srep;

  double as[dim];
  map<key, xyz> geo;
  map<key, double> eff, ori;

  {
    const char * file=getenv("IGEO");
    if(file==NULL) file="../ice/geo-f2k";
    cerr<<"geo file = "<<file<<endl;

    ifstream inFile(file, ifstream::in);
    if(!inFile.fail()){
      string in;
      while(getline(inFile, in)){
	istringstream s(in);
	string mbid, domid;
	int str, dom;
	double x, y, z;
	s>>mbid>>domid>>x>>y>>z>>str>>dom;
	geo[make_pair(str, dom)]=xyz(x, y, z);
      }
      inFile.close();
    }
    cerr<<"geo size = "<<geo.size()<<endl;
  }

  {
    const char * file=getenv("IEFF");
    if(file==NULL) file="../ice/eff-f2k";
    cerr<<"eff file = "<<file<<endl;

    ifstream inFile(file, ifstream::in);
    if(!inFile.fail()){
      string in;
      while(getline(inFile, in)){
	istringstream s(in);
	int str, dom;
	double e;
	s>>str>>dom>>e;
	eff[make_pair(str, dom)]=e;
      }
      inFile.close();
    }
    cerr<<"eff size = "<<eff.size()<<endl;
  }

  {
    const char * file=getenv("IORI");
    if(file==NULL) file="../ice/eff-f2k.ori";
    cerr<<"ori file = "<<file<<endl;

    ifstream inFile(file, ifstream::in);
    if(!inFile.fail()){
      string in;
      while(getline(inFile, in)){
	istringstream s(in);
	int str, dom;
	double e;
	s>>str>>dom>>e;
	ori[make_pair(str, dom)]=e;
      }
      inFile.close();
    }
    cerr<<"ori size = "<<ori.size()<<endl;
  }

  {
    const char * file=getenv("IANG");
    if(file==NULL) file="../ice/as.dat";
    cerr<<"as.dat file = "<<file<<endl;

    int np=-1;
    vector<double> px;

    {
      ifstream inFile(file, ifstream::in);
      if(!inFile.fail()){
	string in;
	while(getline(inFile, in)){
	  istringstream s(in);
	  double x; s>>x;
	  if(np++>=0) px.push_back(x);
	}
	inFile.close();
      }
      cerr<<"as.dat size = "<<np<<endl;
    }

    if(np>0){
      double x=-1, ip=0, y=1, step=2./dim;
      for(int j=0; j<np; j++){ y*=x; ip+=px[j]*y/(j+1); }
      for(int k=0; k<dim; k++){
	double jp=ip;
	x+=step; ip=0, y=1;
	for(int j=0; j<np; j++){ y*=x; ip+=px[j]*y/(j+1); }
	as[k]=(ip-jp)/step;
      }
    }
    else{
      for(int k=0; k<dim; k++) as[k]=1;
    }
  }

  int mode=1;
  if(arg_c>2){
    int aux=atoi(arg_a[2]);
    if(0<=aux && aux<=3) mode=aux;
    else{
      cerr<<"Cannot set mode to "<<aux<<". Keeping the default mode="<<mode<<"."<<endl;
      cerr<<"\t mode=0 \t DOM efficiency"<<endl;
      cerr<<"\t mode=1 \t linear regression"<<endl;
      cerr<<"\t mode=2 \t binned by zenith to emitter"<<endl;
      cerr<<"\t mode=3 \t full likelihood (all terms)"<<endl;
    }
  }

  int num=0;
  vector<double> A, B, C;
  map<key, vector<xyz> > veff;

  if(mode==2){
    num=40;
    A.resize(num*dim), B.resize(num), C.resize(num);
  }

  bool xini=false;
  double x[dim]={0};
  double XX[dim][dim]={0}, XY[dim]={0};

  {
    for(int i=0; i<dim; i++) x[i]=1;
    const char * file=getenv("XINI");
    if(file!=NULL){
      cerr<<"ini file = "<<file<<endl;
      int count=0;

      ifstream inFile(file, ifstream::in);
      if(!inFile.fail()){
	string in;
	while(getline(inFile, in)){
	  istringstream s(in);
	  double c, a, r;
	  s>>c>>a>>r;
	  if(count<dim) x[count++]=r;
	}
	inFile.close();
      }
      cerr<<"ini size = "<<count<<endl;
      xini=true;
    }
  }

  double xmax=0, xsig=0.1;

  {
    char * bmp=getenv("XMAX");
    if(bmp!=NULL){
      xmax=atof(bmp);
      cerr<<"XMAX="<<xmax<<endl;
    }
  }
  {
    char * bmp=getenv("XSIG");
    if(bmp!=NULL){
      xsig=atof(bmp);
      cerr<<"XSIG="<<xsig<<endl;
    }
  }

  bool flag=false;

  if(arg_c>1){
    ifstream inFile(arg_a[1], ifstream::in);
    if(!inFile.fail()){
      int count=0;
      string in;
      while(getline(inFile, in)){
	istringstream s(in);
	int fs, fd, rs, rd;
	double dat, ns, sim[dim], stot=0;
	s>>fs>>fd>>rs>>rd>>ns>>dat;
	for(int i=0; i<dim; i++){ s>>sim[i]; stot+=sim[i]; }

	xyz f=geo[make_pair(fs, fd)], r=geo[make_pair(rs, rd)];
	double dx=f.x-r.x, dy=f.y-r.y, dz=f.z-r.z;
	double dr=sqrt(dx*dx+dy*dy+dz*dz);

	if(stot+dat>0 && dr<350 && rs>0){
	  switch(mode){
	  case 0:
	    {
	      veff[make_pair(rs, rd)].push_back(xyz(stot, dat, ns));
	    }
	    break;

	  case 1:
	    {
	      double err=1;
	      for(int i=0; i<dim; i++) sim[i]/=err;
	      for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++) XX[i][j]+=sim[i]*sim[j];
		XY[i]+=sim[i]*dat/err;
	      }
	    }
	    break;

	  case 2:
	    {
	      int k=min(num-1, max(0, (int) floor(num*(1+dz/dr)/2)));
	      for(int i=0; i<dim; i++) A[k*dim+i]+=sim[i];
	      B[k]+=dat; C[k]+=ns;
	    }
	    break;

	  case 3:
	    {
	      for(int i=0; i<dim; i++) A.push_back(sim[i]);
	      B.push_back(dat); C.push_back(ns); num++;
	    }
	    break;

	  default:
	    ;
	  }

	  count++;
	}
      }
      inFile.close();
      flag=true;
      cerr<<"Accepted entries: "<<count<<endl;
    }
    else{
      cerr<<"Failed on file "<<arg_a[1]<<endl;
    }
  }
  else if(arg_c>0){
    cerr<<"Use: IGEO=[geo-f2k] IEFF=[eff-f2k] IANG=[as.dat] XINI=[x.ini] "<<arg_a[0]<<" [file] [mode]"<<endl;
  }

  if(flag){
    switch(mode){
    case 0:
      for(map<key, vector<xyz> >::const_iterator i=veff.begin(); i!=veff.end(); ++i){
	const key & dom = i->first;
	const vector<xyz> & vset = i->second;
	const int num=vset.size();
	double * a = new double[num];
	double * b = new double[num];
	double * n = new double[num];
	double x[1]={dsr};
	for(int j=0; j<num; j++){
	  a[j]=vset[j].x*srep;
	  b[j]=vset[j].y*drep;
	  n[j]=vset[j].z*noise*drep;
	}

	double e=eff.find(dom)!=eff.end()?eff[dom]:1; e*=escl;
	double o=ori.find(dom)!=ori.end()?ori[dom]:1;
	if(xmax>0){
	  double aux[1]={e}; x[0]*=o/e;
	  double xa[1]={x[0]}, xs[1]={xsig};
	  LSSL::wref(a, b, x, n, num, 1, 0, aux, dsr*o/xmax, dsr*o*xmax, xa, xs);
	}
	else{
	  // for(double y=0.5; y<2.5; y*=1.01){ x[0]=y*dsr; cout<<y<<" "<<LSSL::wllh(a, b, x, n, num, 1)<<endl; }
	  // NNLS::nnls(a, b, x, num, 1);  // x[0]*=o/e;
	  LSSL::wref(a, b, x, n, num, 1);
	}
	e*=x[0]/dsr;
	delete n;
	delete b;
	delete a;
	cout<<dom.first<<" "<<dom.second<<" "<<e<<endl;
      }
    break;

    case 1:
      {
	gsl_matrix * A = gsl_matrix_alloc(dim, dim);
	gsl_vector * B = gsl_vector_alloc(dim);

	for(int i=0; i<dim; i++){
	  for(int j=0; j<dim; j++) gsl_matrix_set(A, i, j, XX[i][j]);
	  gsl_vector_set(B, i, XY[i]);
	}

	gsl_vector * X = gsl_vector_alloc(dim);

	gsl_linalg_HH_solve(A, B, X);

	gsl_matrix_free(A);
	gsl_vector_free(B);

	for(int i=0; i<dim; i++) x[i]=gsl_vector_get(X, i);
	gsl_vector_free(X);
      }
      break;

    case 2:
    case 3:
      {
	double * a = new double[dim*num];
	double * b = new double[num];
	double * n = new double[num];
	for(int i=0; i<dim; i++) x[i]*=dsr;
	for(int j=0; j<num; j++){
	  b[j]=B[j]*drep; n[j]=C[j]*noise*drep;
	  for(int i=0; i<dim; i++) a[i*num+j]=A[i+j*dim]*srep;
	}
	vector<double>().swap(A);
	vector<double>().swap(B);
	vector<double>().swap(C);
	cerr<<"cache cleared"<<endl;

	if(xini){
	  // NNLS::nsum(a, b, x, num, dim);
	  LSSL::wllh(a, b, x, n, num, dim);
	}
	else{
	  if(xmax>0){
	    LSSL::wref(a, b, x, n, num, dim, 0, as, 0, dsr*xmax);
	  }
	  else{
	    // NNLS::nnls(a, b, x, num, dim);
	    LSSL::wref(a, b, x, n, num, dim);
	  }
	}
	for(int i=0; i<dim; i++) x[i]/=dsr;
	delete n;
	delete b;
	delete a;
      }
      break;

    default:
      break;
    }

    if(mode>0 && !xini){
      for(int i=0; i<dim; i++) cout<<((2./dim)*(i+0.5)-1)<<" "<<as[i]<<" "<<x[i]<<endl;
    }
  }
}
