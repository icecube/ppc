namespace LSSL{
  const bool verbose=true;
  const double xx=1.e-8;
  double SIG=0.10;

  int n, m;
  double *a, *c, *d, *s;

  double llh_stat(double *x, double *dx){
    const int iter=100;
    double sum=0;
    if(dx!=NULL) for(int l=0; l<m; l++) dx[l]=0;

    for(int k=0; k<n; k++){
      double dk=d[k], eta=s[k];

      double xi=0, xm=0;

      for(int l=0; l<m; l++){
	double wl=x[l], s=a[k+n*l];
	if(s>0 && wl>0){
	  double xv=-1/wl;
	  if(xm==0 || xv>xm) xm=xv;
	}
      }

      if(dk==0) xi=1;
      else if(xm==0) xi=1-dk/eta;
      else{
	double xo=xi;

	for(int i=0; i<iter; i++){
	  double f=eta, df=0;
	  for(int l=0; l<m; l++){
	    double wl=x[l], s=a[k+n*l];
	    if(s>0){
	      double r=wl/(1+xi*wl);
	      f+=s*r; df+=s*r*r;
	    }
	  }
	  df=df/(f*f)+1./dk;
	  f=1/f-(1-xi)/dk;

	  {
	    double xn=xi-f/df;
	    if(xn<=xm) xn=(xi+xm)/2;
	    else if(xn>=1) xn=(1+xi)/2;
	    xi=xn;
	  }

	  if(fabs(xi-xo)<=fabs(xi+xo)*xx) break;
	  xo=xi;
	}
      }

      double sk=0;
      for(int l=0; l<m; l++){
	double wl=x[l], s=a[k+n*l];
	if(s>0){
	  double aux=max(xx, 1+xi*wl); sk+=s*log(aux);
	  if(dx!=NULL) dx[l]+=s*xi/aux;
	}
      }
      if(xi<1) sk+=dk*log(1-xi);
      sum+=sk;
    }

    return sum;
  }

  double llh_syst(double *x, double *dx){
    const int iter=100;
    const double s2=SIG*SIG;
    double sum=0;
    if(dx!=NULL) for(int l=0; l<m; l++) dx[l]=0;

    for(int k=0; k<n; k++){
      double dk=d[k], eta=s[k];
      double ks=dk*s2, es=log(eta*s2);
      double wx=0, wn=-ks, zi;

      for(int l=0; l<m; l++){
	double wl=x[l], s=a[k+n*l];
	double sw=s*wl;
	if(sw>0) if(wx==0 || wl>wx) wx=wl;
      }

      {
	double z=es;
	for(int i=0; i<iter; i++){
	  double zk=z+ks;
	  double f=zk>1?z+log(zk)-es:zk-exp(es-z);
	  double df=zk>1?1+1/zk:1+zk;

	  double dz=f/df;
	  double zn=z-dz;
	  if(zn>wn) z=zn; else z=(z+wn)/2;
	  if(fabs(dz)<xx) break;
	}
	zi=z;
      }

      if(wx>0){
	double z=0;
	for(int i=0; i<iter; i++){
	  double exz=exp(-z);
	  double f=(z+ks)+z*wx*exz;
	  double df=(1+z+ks)+wx*exz;

	  double dz=f/df;
	  double zn=z-dz;
	  if(zn>wn) z=zn; else z=(z+wn)/2;
	  if(fabs(dz)<xx) break;
	}

	wx=max(z, zi);
	zi=max(zi+1, 0.);

	for(int i=0; i<iter; i++){
	  double exz=exp(zi);
	  double mu=exz*(dk+zi/s2);
	  double mp=exz*(dk+(1+zi)/s2);

	  double f=1-eta/mu, df=eta*mp/(mu*mu);

	  for(int l=0; l<m; l++){
	    double wl=x[l], s=a[k+n*l];
	    double sw=s*wl;
	    if(sw>0){
	      double r=mu+wl*zi/s2;
	      f-=sw/r; df+=(mp+wl/s2)*sw/(r*r);
	    }
	  }
	  double dz=f/df;

	  {
	    double zn=zi-dz;
	    if(zn<=wx) zn=(zi+wx)/2;
	    zi=zn;
	  }

	  // cerr<<i<<" "<<zi<<" "<<dz<<endl;
	  if(fabs(dz)<xx) break;
	}
      }

      {
	double sk=0;
	double exz=exp(zi);
	double md=dk+zi/s2;
	double mu=exz*md;

	for(int l=0; l<m; l++){
	  double wl=x[l], s=a[k+n*l];
	  double sw=s*wl;
	  double di=zi/(s2*mu);

	  if(sw>0){
	    double aux=max(xx, 1+wl*di); sk+=s*log(aux);
	    if(dx!=NULL) dx[l]+=s*di/aux;
	  }
	}
	if(dk>0) sk+=dk*log(dk/md);
	sk+=zi*zi/(2*s2);
	sum+=sk;
      }
    }

    return sum;
  }

  double *xa, *xs;

  double llh(double *x, double *dx = NULL){  // negative log of likelihood ratio
    double sum=SIG>0?llh_syst(x, dx):llh_stat(x, dx);
    if(!isfinite(sum)){
      cerr<<"Error: llh = nan or inf!"<<endl;
      sum=-n*log(xx);
    }

    if(xs!=NULL){
      double sub=0;
      for(int i=0; i<m; i++){
	double aux=log(x[i]/xa[i])/xs[i];
	sub+=aux*aux/2;
	if(dx!=NULL) dx[i]+=aux/(x[i]*xs[i]);
      }
      sum+=sub;
    }
    /*
    else if(xa!=NULL){
      double sub=0, cs=1.e-3;
      for(int i=0; i<m; i++){
	sub+=xa[i]*x[i];
	if(dx!=NULL) if(x[i]>0) dx[i]+=cs*xa[i];
      }
      sum+=cs*sub;
    }
    */
    return sum;
  }

  double lnllh(double *x, double *dx = NULL){  // negative log of likelihood ratio
    double X[m];
    for(int i=0; i<m; i++) X[i]=exp(x[i]);
    double f=llh(X, dx);
    if(dx!=NULL) for(int i=0; i<m; i++) dx[i]*=X[i];
    return f;
  }

  void wref_fdf(const gsl_vector *x, void *par, double *f, gsl_vector *df){
    double X[m], D[m];
    for(int i=0; i<m; i++) X[i]=exp(gsl_vector_get(x, i))/c[i];
    double result=llh(X, D); // cerr<<result<<endl;

    if(f!=NULL) *f = result;
    if(df!=NULL) for(int i=0; i<m; i++) gsl_vector_set(df, i, X[i]*D[i]);
  }

  void wref_df(const gsl_vector *x, void *par, gsl_vector *df){
    wref_fdf(x, par, NULL, df);
  }

  double wref_f(const gsl_vector *x, void *par){
    double result;
    wref_fdf(x, par, &result, NULL);
    return result;
  }

  double wllh(double * A, double * D, double * X, double * S, int N, int M){
    n=N; m=M; a=A; d=D; s=S;
    double F=llh(X); if(verbose) cerr<<"X0 "<<F<<" *"<<endl;
    return F/N;
  }

  double wref(double * A, double * D, double * V, double * S, int N, int M, double trs = 0, double *as = NULL, double xmin = 0, double xmax = 1, double *XA = NULL, double *XS = NULL){
    n=N; m=0;
    bool R[M];
    double stot=0, dtot=0;
    for(int i=0; i<M; i++){
      double ci=0;
      for(int j=0; j<N; j++) ci+=A[N*i+j];
      if(trs<-1) V[i]=0;
      if(R[i] = ci>0 && V[i]>trs) m++;
      stot+=ci;
    }

    if(m==0) return llh(NULL)/N;
    for(int j=0; j<N; j++) dtot+=D[j];
    double Xfl=dtot/stot;

    double X[m]; a = new double[n*m]; c = new double[m]; d=D; s=S;
    if(XA!=NULL) xa = new double[m]; else xa = NULL;
    if(XS!=NULL) xs = new double[m]; else xs = NULL;

    for(int i=0, k=0; i<M; i++){
      if(R[i]){
	for(int j=0; j<N; j++) a[n*k+j]=A[N*i+j];
	if(XA!=NULL) xa[k]=XA[i];
	if(XS!=NULL) xs[k]=XS[i];
	X[k++]=V[i];
      }
    }

    for(int i=0; i<m; i++){
      double ci=0;
      for(int j=0; j<n; j++) ci+=a[j+n*i];
      c[i]=ci;
      if(trs<-1 || trs<0&&X[i]<=0) X[i]=Xfl;
    }

    // n=N, m=M; a=A, d=D;  // double A[(0...n)+n*(0...m)], D[n], X[m];
    double result=NAN;

    if(as!=NULL){
      { // NMML
	double x[2][m], g[2][m];
	for(int i=0; i<m; i++) x[0][i]=X[i];

	double F=llh(x[0], g[0]); if(verbose) cerr<<"X0 "<<F<<" *"<<endl;
	result=F;

	for(int i=0; i<m; i++){
	  const double xl=xmin/as[i];
	  const double xh=xmax/as[i];
	  if(x[0][i]<=xl && g[0][i]>0){
	    x[0][i]=xl; g[0][i]=0;
	  }
	  else if(x[0][i]>=xh && g[0][i]<0){
	    x[0][i]=xh; g[0][i]=0;
	  }
	  else{
	    g[0][i]*=x[0][i]/c[i];
	  }
	}

	double alf=1;
	int k0, k1;
	for(int k=1; k<100; k++){
	  k0=(k-1)%2, k1=k%2;
	  bool flag=true;
	  double F1;
	  while(flag){
	    for(int i=0; i<m; i++){
	      const double xl=xmin/as[i];
	      const double xh=xmax/as[i];
	      x[k1][i]=min(xh, max(xl, x[k0][i]-alf*g[k0][i]));
	    }
	    F1=llh(x[k1], g[k1]);
	    flag = F1>F && alf>1;
	    if(flag) alf/=2; if(alf<1) alf=1;
	  }
	  F=F1;

	  double dot=0; alf=0;
	  for(int i=0; i<m; i++){
	    double xi=0, di=0;
	    const double xl=xmin/as[i];
	    const double xh=xmax/as[i];
	    if(x[k1][i]<=xl && g[k1][i]>0){
	      x[k1][i]=xl; g[k1][i]=0;
	    }
	    else if(x[k1][i]>=xh && g[k1][i]<0){
	      x[k1][i]=xh; g[k1][i]=0;
	    }
	    else{
	      g[k1][i]*=x[k1][i]/c[i];
	      xi=x[k1][i]-x[k0][i];
	      di=g[k1][i]-g[k0][i];
	      alf+=xi*xi, dot+=xi*di;
	    }
	  }
	  alf/=dot;
	  if(!isfinite(alf) || alf<1) alf=1;
	}

	int use=F>result?0:F==result?1:2; // use=3;
	if(use>0){
	  result=F;
	  for(int i=0; i<m; i++) X[i]=x[k1][i];
	}
	if(verbose) cerr<<"XF"<<" "<<F<<" "<<(use>1?"*":use>0?"+":"-")<<endl;
      }
    }

    else{
      double x=SIG; SIG=0;
      for(int k=0, n=x>0?2:1; k<n; k++, SIG=x){ // NMML
	double x[2][m], g[2][m];
	for(int i=0; i<m; i++) x[0][i]=X[i];

	double F=llh(x[0], g[0]); if(verbose) cerr<<"X0 "<<F<<" *"<<endl;
	result=F;

	for(int i=0; i<m; i++){
	  if(x[0][i]==0 && g[0][i]>0){
	    g[0][i]=0;
	  }
	  else{
	    g[0][i]*=x[0][i]/c[i];
	  }
	}

	double alf=1;
	int k0, k1;
	for(int k=1; k<100; k++){
	  k0=(k-1)%2, k1=k%2;
	  bool flag=true;
	  double F1;
	  while(flag){
	    for(int i=0; i<m; i++) x[k1][i]=max(0., x[k0][i]-alf*g[k0][i]);
	    F1=llh(x[k1], g[k1]);
	    flag = F1>F && alf>1;
	    if(flag) alf/=2; if(alf<1) alf=1;
	  }
	  F=F1;

	  double dot=0; alf=0;
	  for(int i=0; i<m; i++){
	    double xi=0, di=0;
	    if(x[k1][i]==0 && g[k1][i]>0){
	      g[k1][i]=0;
	    }
	    else{
	      g[k1][i]*=x[k1][i]/c[i];
	      xi=x[k1][i]-x[k0][i];
	      di=g[k1][i]-g[k0][i];
	      alf+=xi*xi, dot+=xi*di;
	    }
	  }
	  alf/=dot;
	  if(!isfinite(alf) || alf<1) alf=1;
	}

	int use=F>result?0:F==result?1:2;
	if(use>0){
	  result=F;
	  for(int i=0; i<m; i++) X[i]=x[k1][i];
	}
	if(verbose) cerr<<"XF"<<" "<<F<<" "<<(use>1?"*":use>0?"+":"-")<<endl;
      }

      for(int k=0; k<4; k++){
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf G;

	G.n = m;
	G.f = wref_f;
	G.df = wref_df;
	G.fdf = wref_fdf;
	G.params = NULL;

	x = gsl_vector_alloc(m);
	for(int i=0; i<m; i++) gsl_vector_set(x, i, log(X[i]*c[i]));

	switch(k){
	case 0: T = gsl_multimin_fdfminimizer_conjugate_fr; break;
	case 1: T = gsl_multimin_fdfminimizer_conjugate_pr; break;
	case 2: T = gsl_multimin_fdfminimizer_vector_bfgs2; break;
	default: T = gsl_multimin_fdfminimizer_steepest_descent;
	}
	s = gsl_multimin_fdfminimizer_alloc(T, m);

	gsl_multimin_fdfminimizer_set(s, &G, x, 0.01, 0.1);

	for(int iter=0; iter<25; iter++) if(gsl_multimin_fdfminimizer_iterate(s)!=GSL_SUCCESS) break;

	double F=s->f;

	int use=F>result?0:F==result?1:2;
	if(use>0){
	  result=F;
	  for(int i=0; i<m; i++) X[i]=exp(gsl_vector_get(s->x, i))/c[i];
	}
	if(verbose) cerr<<"Z"<<k<<" "<<F<<" "<<(use>1?"*":use>0?"+":"-")<<endl;

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);
      }
    }

    for(int i=0, k=0; i<M; i++) if(R[i]) V[i]=X[k++];

    if(XA!=NULL) delete xa;
    if(XS!=NULL) delete xs;

    delete a; delete c;
    return result/N;
  }
}
