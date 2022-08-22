#include <map>
#include <set>
#include <vector>

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

struct dir{
  double n[3];

  dir(){}

  dir(const double n[]){
    for(int i=0; i<3; i++) this->n[i]=n[i];
  }

  bool operator< (const dir& rhs) const {
    return n[0] < rhs.n[0] || (n[0] == rhs.n[0] && (n[1] < rhs.n[1] || (n[1] == rhs.n[1] && n[2] < rhs.n[2])));
  }
};

struct vert:dir{
  double f;

  vert(){}

  vert(const dir n, double f) : dir(n){
    this->f=f;
  }
};

ostream& operator<<(ostream& out, const dir& n){
  out<<n.n[0]<<" "<<n.n[1]<<" "<<n.n[2];
  return out;
}

struct face{
  dir v[3];

  face(const dir & x, const dir & y, const dir & z){
    v[0]=x, v[1]=y, v[2]=z;
  }

  void half(dir h[]){
    for(int i=0; i<3; i++){
      h[0].n[i]=v[1].n[i]+v[2].n[i];
      h[1].n[i]=v[2].n[i]+v[0].n[i];
      h[2].n[i]=v[0].n[i]+v[1].n[i];
    }
    double r[3]={0, 0, 0};
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) r[i]+=h[i].n[j]*h[i].n[j];
    for(int i=0; i<3; i++){
      r[i]=1/sqrt(r[i]);
      for(int j=0; j<3; j++) h[i].n[j]*=r[i];
    }
  }
};

struct Ico{
  int vn, fn;
  vector<dir> v;
  vector< vector<int> > f;

  Ico(int num) : vn(0), fn(0){
    vector<face> set;

    {
      const double t=(1+sqrt(5))/2;
      const double r=sqrt(1+t*t);
      const double a=t/r, b=1/r;

      const int vn=12, fn=20;
      double v[vn][3];
      int f[fn][3];

      for(int i=0; i<vn; i++){
	int x=i/4, j=i%4;
	int y=(x+1)%3, z=(x+2)%3;
	v[i][x]=j==0||j==3?a:-a;
	v[i][y]=j==0||j==1?b:-b;
	v[i][z]=0;
      }

      for(int i=0, n=0; i<vn; i++)
	for(int j=i+1; j<11; j++)
	  if(v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2]>0)
	    for(int k=j+1; k<vn; k++)
	      if(v[j][0]*v[k][0]+v[j][1]*v[k][1]+v[j][2]*v[k][2]>0 &&
		 v[k][0]*v[i][0]+v[k][1]*v[i][1]+v[k][2]*v[i][2]>0) f[n][0]=i, f[n][1]=j, f[n++][2]=k;

      if(false){
	for(int i=0; i<vn; i++) cout<<v[i][0]<<" "<<v[i][1]<<" "<<v[i][2]<<endl;
	for(int i=0; i<fn; i++) cout<<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<endl;
      }

      for(int i=0; i<fn; i++){
	int x=f[i][0], y=f[i][1], z=f[i][2];
	set.push_back(face(dir(v[x]), dir(v[y]), dir(v[z])));
      }
    }

    for(int k=0; k<num; k++){
      vector<face> add;

      for(vector<face>::iterator i=set.begin(); i!=set.end(); i++){
	face & one = *i;
	dir * w = one.v;

	dir u[3];
	one.half(u);

	add.push_back(face(w[0], u[1], u[2]));
	add.push_back(face(w[1], u[2], u[0]));
	add.push_back(face(w[2], u[0], u[1]));
	add.push_back(face(u[0], u[1], u[2]));
      }

      set=add;
    }

    {
      map<dir, int> dirs;
      for(vector<face>::const_iterator i=set.begin(); i!=set.end(); i++){
	const dir * w = i->v;
	for(int i=0; i<3; i++){
	  const dir & wi = w[i];
	  if(dirs.find(wi)==dirs.end()) dirs.insert(make_pair(wi, 0));
	}
      }

      for(map<dir, int>::iterator i=dirs.begin(); i!=dirs.end(); i++) v.push_back(i->first), i->second=vn++;
      for(vector<face>::const_iterator i=set.begin(); i!=set.end(); i++){
	vector<int> tmp;
	for(int j=0; j<3; j++) tmp.push_back(dirs[i->v[j]]);
	f.push_back(tmp); fn++;
      }
    }

    cerr<<"Icosahedral mesh("<<num<<"): points="<<vn<<", triangles="<<fn<<endl;
  }
};

main(int arg_c, char *arg_a[]){
  int n=-1;
  if(arg_c>1) n=atoi(arg_a[1]);
  Ico ico(n);
  for(int i=0; i<ico.vn; i++) cout<<i<<" "<<ico.v[i]<<endl;
}
