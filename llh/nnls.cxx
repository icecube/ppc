namespace NNLS{

  /*
    Non-negative least squares by C. Lawson and R. Hanson.
  */

  double sign(double x, double y){
    return y>=0?fabs(x):-fabs(x);
  }

  double square(double x){
    return x*x;
  }

  void g1(double a, double b, double & cterm, double & sterm, double & sig){
    const double one=1, zero=0;
    double xr, yr;

    if(fabs(a)>fabs(b)){
      xr=b/a;
      yr=sqrt(one+xr*xr);
      cterm=sign(one/yr,a);
      sterm=cterm*xr;
      sig=fabs(a)*yr;
    }
    else if(b!=zero){
      xr=a/b;
      yr=sqrt(one+xr*xr);
      sterm=sign(one/yr,b);
      cterm=sterm*xr;
      sig=fabs(b)*yr;
    }
    else{
      sig=zero;
      cterm=zero;
      sterm=one;
    }
  }

  void h12(int mode, int lpivot, int l1, int m, double u[], int iue,
	   double & up, double c[], int ice, int icv, int ncv){
    const double one=1;
    int i, i2, i3, i4, incr, j;
    double b, cl, clinv, sm;

#define c(x) c[x-1]
#define u(x,y) u[x-1+(y-1)*iue]

    if(0>=lpivot || lpivot>=l1 || l1>m) return;
    cl=fabs(u(1,lpivot));

    if(mode==1){
      for(j=l1; j<=m; j++) cl=fmax(fabs(u(1,j)),cl);
      if (cl<=0) return;
      clinv=one/cl;
      sm=square(u(1,lpivot)*clinv);
      for(j=l1; j<=m; j++) sm+=square(u(1,j)*clinv);
      cl*=sqrt(sm);
      if(u(1,lpivot)>0) cl=-cl;
      up=u(1,lpivot)-cl;
      u(1,lpivot)=cl;
    }
    else if(cl<=0) return;
    if(ncv<=0) return;
    b=up*u(1,lpivot);
    if(b>=0) return;

    b=one/b;
    i2=1-icv+ice*(lpivot-1);
    incr=ice*(l1-lpivot);

    for(j=1; j<=ncv; j++){
      i2+=icv;
      i3=i2+incr;
      i4=i3;
      sm=c(i2)*up;
      for(i=l1; i<=m; i++){
	sm+=c(i3)*u(1,i);
	i3+=ice;
      }
      if(sm!=0){
	sm*=b;
	c(i2)+=sm*up;
	for(i=l1; i<=m; i++){
	  c(i4)+=sm*u(1,i);
	  i4+=ice;
	}
      }
    }

#undef c
#undef u

    return;
  }

  int nnls(double a[], int mda, int m, int n, double b[], double x[],
	   double & rnorm, double w[], double zz[], int index[], int & nsetp){
    int i, ii, ip, iter, itmax, iz, iz1, iz2, izmax, j, jj, jz, l, npp1, mode;
    double alpha, asave, cc, dummy, sm, ss, t, temp, unorm, up, wmax, ztest;
    const double zero=0, two=2, factor=0.01;

#define b(x) b[x-1]
#define x(y) x[y-1]
#define w(x) w[x-1]
#define zz(x) zz[x-1]
#define index(x) index[x-1]
#define a(x,y) a[x-1+(y-1)*mda]

    if(m<=0 || n<=0) return 2;
    mode=1;
    iter=0;
    itmax=3*n;

    for(i=1; i<=n; i++){
      x(i)=zero;
      index(i)=i;
    }

    iz2=n;
    iz1=1;
    nsetp=0;
    npp1=1;

    while(true){
      if(iz1>iz2 || nsetp>=m) break;

      for(iz=iz1; iz<=iz2; iz++){
	j=index(iz);
	sm=zero;
	for(l=npp1; l<=m; l++) sm+=a(l,j)*b(l);
	w(j)=sm;
      }

      bool flag=false;
      while(true){
	wmax=zero;
	for(iz=iz1; iz<=iz2; iz++){
	  j=index(iz);
	  if(w(j)>wmax){
	    wmax=w(j);
	    izmax=iz;
	  }
	}

	if(wmax<=zero){ flag=true; break; }
	iz=izmax;
	j=index(iz);

	asave=a(npp1,j);
	h12(1, npp1, npp1+1, m, &a(1,j), 1, up, &dummy, 1, 1, 0);
	unorm=zero;
	if(nsetp!=0){
	  for(l=1; l<=nsetp; l++) unorm+=square(a(l,j));
	}
	unorm=sqrt(unorm);

	if(unorm+fabs(a(npp1,j))*factor>unorm){
	  for(l=1; l<=m; l++) zz(l)=b(l);
	  h12(2, npp1, npp1+1, m, &a(1,j), 1, up, zz, 1, 1, 1);
	  ztest=zz(npp1)/a(npp1,j);
	  if(ztest>zero) break;
	}

	a(npp1,j)=asave;
	w(j)=zero;
      }
      if(flag) break;

      for(l=1; l<=m; l++) b(l)=zz(l);

      index(iz)=index(iz1);
      index(iz1++)=j;
      nsetp=npp1++;

      if(iz1<=iz2){
	for(jz=iz1; jz<=iz2; jz++){
	  jj=index(jz);
	  h12(2, nsetp, npp1, m, &a(1,j), 1, up, &a(1,jj), 1, mda, 1);
	}
      }

      if(nsetp!=m){
	for(l=npp1; l<=m; l++) a(l,j)=zero;
      }
      w(j)=zero;

      for(l=1; l<=nsetp; l++){
	ip=nsetp+1-l;
	if(l!=1){
	  for(ii=1; ii<=ip; ii++) zz(ii)-=a(ii,jj)*zz(ip+1);
	}
	jj=index(ip);
	zz(ip)/=a(ip,jj);
      }

      flag=false;
      while(true){
	// for(int i=0; i<n; i++){ for(int j=0; j<m; j++) cout<<" "<<a(i,j); cout<<endl; } cout<<endl;
	iter++;
	if(iter>itmax){
	  mode=3;
	  cerr<<"nnls quitting on iteration count!"<<endl;
	  flag=true; break;
	}

	alpha=two;
	for(ip=1; ip<=nsetp; ip++){
	  l=index(ip);
	  if(zz(ip)<=zero){
	    t=-x(l)/(zz(ip)-x(l));
	    if(alpha>t){
	      alpha=t;
	      jj=ip;
	    }
	  }
	}

	if(alpha==two) break;

	for(ip=1; ip<=nsetp; ip++){
	  l=index(ip);
	  x(l)+=alpha*(zz(ip)-x(l));
	}

	i=index(jj);

	while(true){
	  x(i)=zero;

	  if(jj!=nsetp){
	    jj++;
	    for(j=jj; j<=nsetp; j++){
	      ii=index(j);
	      index(j-1)=ii;
	      g1(a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii));
	      a(j,ii)=zero;
	      for(l=1; l<=n; l++){
		if(l!=ii){
		  temp=a(j-1,l);
		  a(j-1,l)=cc*temp+ss*a(j,l);
		  a(j,l)=-ss*temp+cc*a(j,l);
		}
	      }
	      temp=b(j-1);
	      b(j-1)=cc*temp+ss*b(j);
	      b(j)=-ss*temp+cc*b(j);
	    }
	  }

	  npp1=nsetp--;
	  index(--iz1)=i;

	  for(jj=1; jj<=nsetp; jj++){
	    i=index(jj);
	    if(x(i)<=zero) continue;
	  }

	  for(i=1; i<=m; i++) zz(i)=b(i);

	  for(l=1; l<=nsetp; l++){
	    ip=nsetp+1-l;
	    if(l!=1){
	      for(ii=1; ii<=ip; ii++) zz(ii)-=a(ii,jj)*zz(ip+1);
	    }
	    jj=index(ip);
	    zz(ip)/=a(ip,jj);
	  }
	  break;
	}
      }
      if(flag) break;

      for(ip=1; ip<=nsetp; ip++){
	i=index(ip);
	x(i)=zz(ip);
      }
    }

    sm=zero;
    if(npp1<=m){
      for(i=npp1; i<=m; i++) sm+=square(b(i));
    }
    else{
      for(j=1; j<=n; j++) w(j)=zero;
    }

    rnorm=sqrt(sm);

#undef b
#undef x
#undef w
#undef zz
#undef index
#undef z

    return mode;
  }

  double nnls(double a[], double b[], double x[], int m, int n){
    int mda=m, nsetp, index[n];
    double rnorm, w[n], zz[m];

    double * B = new double[m];
    double * A = new double[n*m];

    for(int i=0; i<m; i++) B[i]=b[i];
    for(int i=0; i<n*m; i++) A[i]=a[i];

    nnls(A, mda, m, n, B, x, rnorm, w, zz, index, nsetp);

    delete B;
    delete A;
    return rnorm/m;
  }

  double nsum(double a[], double b[], double x[], int m, int n){
    double sum=0;
    for(int i=0; i<m; i++){
      double qi=0;
      for(int j=0; j<n; j++) qi+=a[i+j*m]*x[j];
      double aux=b[i]-qi;
      sum+=aux*aux;
    }
    return sqrt(sum)/m;
  }
}
