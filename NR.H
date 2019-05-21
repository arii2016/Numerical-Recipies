	typedef struct FCOMPLEX {float r,i;} fcomplex;
	typedef struct IMMENSE {unsigned long l,r;} immense;
	typedef struct GREAT {unsigned short l,c,r;} great;

	void  adi(double **a, double **b, double **c, double **d, double **e,
					double **f, double **g, double **u, int jmax, int k,
					double alpha, double beta, double eps);
	void  amoeba(float **p, float *y, int ndim, float ftol,
					float (*funk)(float *), int *iter);
	void  anneal(float *x, float *y, int *iorder, int ncity);
	void  avevar(float *data, int n, float *ave, float *svar);
	void  balanc(float **a, int n);
	void  bcucof(float *y, float *y1, float *y2, float *y12, float d1,
					float d2, float **c);
	void  bcuint(float *y, float *y1, float *y2, float *y12, float x1l,
					float x1u, float x2l, float x2u, float x1, float x2,
					float *ansy, float *ansy1, float *ansy2);
	float bessi(int n, float x);
	float bessi0(float x);
	float bessi1(float x);
	float bessj(int n, float x);
	float bessj0(float x);
	float bessj1(float x);
	float bessk(int n, float x);
	float bessk0(float x);
	float bessk1(float x);
	float bessy(int n, float x);
	float bessy0(float x);
	float bessy1(float x);
	float beta(float z, float w);
	float betacf(float a, float b, float x);
	float betai(float a, float b, float x);
	float bico(int n, int k);
	void  bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
	float bnldev(float pp, int n, int *idum);
	float brent(float ax, float bx, float cx, float (*f)(float), float tol,
					float *xmin);
	void  bsstep(float *y, float *dydx, int nv, float *xx, float htry,
					float eps, float *yscal, float *hdid, float *hnext,
					void (*derivs)(float,float *,float *));
	void  caldat(long julian, int *mm, int *id, int *iyyy);
	float cel(float qqc, float pp, float aa, float bb);
	void  chder(float a, float b, float *c, float *cder, int n);
	float chebev(float a, float b, float *c, int m, float x);
	void  chebft(float a, float b, float *c, int n, float (*func)(float));
	void  chebpc(float *c, float *d, int n);
	void  chint(float a, float b, float *c, float *cint, int n);
	void  chsone(float *bins, float *ebins, int nbins, int knstrn,
					float *df, float *chsq, float *prob);
	void  chstwo(float *bins1, float *bins2, int nbins, int knstrn,
					float *df, float *chsq, float *prob);
	void  cntab1(int **nn, int n1, int nj, float *chisq, float *df,
					float *prob, float *cramrv, float *ccc);
	void  cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy,
					float *hygx, float *hxgy, float *uygx, float *uxgy,
					float *uxy);
	void  convlv(float *data, int n, float *respns, int m, int isign,
					float *ans);
	void  correl(float *data1, float *data2, int n, float *ans);
	void  cosft(float *y, int n, int isign);
	void  covsrt(float **covar, int ma, int *lista, int mfit);
	void  crank(int n, float *w, float *s);
	float dbrent(float ax, float bx, float cx, float (*f)(float),
					float (*df)(float), float tol, float *xmin);
	void  ddpoly(float *c, int nc, float x, float *pd, int nd);
	void  ks(immense key, int n, great *kn);
	void  cyfun(unsigned long ir, great k, unsigned long *iout);
	float df1dim(float x);
	void  dfpmin(float *p, int n, float ftol, int *iter, float *fret,
					float (*func)(float *), void (*dfunc)(float *,float *));
	void  difeq(int k, int k1, int k2, int jsf, int is1, int isf,
					int *indexv, int ne, float **s, float **y);
	void  dlinmin(float *p, float *xi, int n, float *fret,
					float (*func)(float *), void (*dfunc)(float *,float *));
	void  eclass(int *nf, int n, int *lista, int *listb, int m);
	void  eclazz(int *nf, int n, int (*equiv)(int,int));
	void  eigsrt(float *d, float **v, int n);
	float el2(float x, float qqc, float aa, float bb);
	void  elmhes(float **a, int n);
	float erf(float x);
	float erfc(float x);
	float erfcc(float x);
	void  eulsum(float *sum, float term, int jterm, float *wksp);
	float evlmem(float fdt, float *cof, int m, float pm);
	float expdev(int *idum);
	float f1dim(float x);
	float factln(int n);
	float factrl(int n);
	void  fgauss(float x, float *a, float *y, float *dyda, int na);
	void  fit(float *x, float *y, int ndata, float *sig, int mwt, float *a,
					float *b, float *siga, float *sigb, float *chi2, float *q);
	void  fixrts(float *d, int npoles);
	void  fleg(float x, float *pl, int nl);
	void  flmoon(int n, int nph, long *jd, float *frac);
	void  four1(float *data, int nn, int isign);
	void  fourn(float *data, int *nn, int ndim, int isign);
	void  fpoly(float x, float *p, int np);
	void  frprmn(float *p, int n, float ftol, int *iter, float *fret,
					float (*func)(float *), void (*dfunc)(float *,float *));
	void  ftest(float *data1, int n1, float *data2, int n2, float *f,
					float *prob);
	float gamdev(int ia, int *idum);
	float gammln(float xx);
	float gammp(float a, float x);
	float gammq(float a, float x);
	float gasdev(int *idum);
	void  gauleg(double x1, double x2, double *x, double *w, int n);
	void  gaussj(float **a, int n, float **b, int m);
	void  gcf(float *gammcf, float a, float x, float *gln);
	float golden(float ax, float bx, float cx, float (*f)(float), float tol,
					float *xmin);
	void  gser(float *gamser, float a, float x, float *gln);
	void  hqr(float **a, int n, float *wr, float *wi);
	void  hunt(float *xx, int n, float x, int *jlo);
	void  indexx(int n, float *arrin, int *indx);
	int   irbit1(unsigned long int *iseed);
	int   irbit2(unsigned long int *iseed);
	void  jacobi(float **a, int n, float *d, float **v, int *nrot);
	long  julday(int mm, int id, int iyyy);
	void  kendl1(float *data1, float *data2, int n, float *tau, float *z,
					float *prob);
	void  kendl2(float **tab, int i, int j, float *tau, float *z,
					float *prob);
	void  ksone(float *data, int n, float (*func)(float), float *d,
					float *prob);
	void  kstwo(float *data1, int n1, float *data2, int n2, float *d,
					float *prob);
	void  laguer(fcomplex *a, int m, fcomplex *x, float eps, int polish);
	void  lfit(float *x, float *y, float *sig, int ndata, float *a, int ma,
					int *lista, int mfit, float **covar, float *chisq,
					void (*funcs)(float,float *,int));
	void  linmin(float *p, float *xi, int n, float *fret, float (*func)(float));
	void  locate(float *xx, int n, float x, int *j);
	void  lubksb(float **a, int n, int *indx, float *b);
	void  ludcmp(float **a, int n, int *indx, float *d);
	void  mdian1(float *x, int n, float *xmed);
	void  mdian2(float *x, int n, float *xmed);
	void  medfit(float *x, float *y, int ndata, float *a, float *b,
					float *abdev);
	void  memcof(float *data, int n, int m, float *pm, float *cof);
	float midexp(float (*funk)(float), float aa, float bb, int n);
	float midinf(float (*funk)(float), float aa, float bb, int n);
	float midpnt(float (*func)(float), float a, float b, int n);
	float midsql(float (*funk)(float), float aa, float bb, int n);
	float midsqu(float (*funk)(float), float aa, float bb, int n);
	void  mmid(float *y, float *dydx, int nvar, float xs, float htot,
					int nstep, float *yout,
					void (*derivs)(float,float *,float *));
	void  mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
					float *fc, float (*func)(float));
	void  mnewt(int ntrial, float *x, int n, float tolx, float tolf);
	void  moment(float *data, int n, float *ave, float *adev, float *sdev,
					float *svar, float *skew, float *curt);
	void  mprove(float **a, float **alud, int n, int *indx, float *b,
					float *x);
	void  mrqcof(float *x, float *y, float *sig, int ndata, float *a, int ma,
					int *lista, int mfit, float **alpha, float *beta, float
					*chisq, void (*funcs)(float,float *,float *,float *,int));
	void  mrqmin(float *x, float *y, float *sig, int ndata, float *a,
					int ma, int *lista, int mfit, float **covar, float **alpha,
					float *chisq, void (*funcs)(float,float *,float *,float *,int),
					float *alamda);
	void  odeint(float *ystart, int nvar, float x1, float x2, float eps,
					float h1, float hmin, int *nok, int *nbad,
					void (*derivs)(float,float *,float *),
					void  (*rkqc)(float *,float *,int,float *,float,float,float *,
					float *,float *,void (*)(float,float *,float *)));
	void  pcshft(float a, float b, float *d, int n);
	void  pearsn(float *x, float *y, int n, float *r, float *prob, float *z);
	void  piksr2(int n, float *arr, float *brr);
	void  piksrt(int n, float *arr);
	void  pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
					float ***c, float **s);
	float plgndr(int l, int m, float x);
	float poidev(float xm, int *idum);
	void  polcoe(float *x, float *y, int n, float *cof);
	void  polcof(float *xa, float *ya, int n, float *cof);
	void  poldiv(float *u, int n, float *v, int nv, float *q, float *r);
	void  polin2(float *x1a, float *x2a, float **ya, int m, int n, float x1,
					float x2, float *y, float *dy);
	void  polint(float *xa, float *ya, int n, float x, float *y, float *dy);
	void  powell(float *p, float **xi, int n, float ftol, int *iter,
					float *fret, float (*func)(float *));
	void  predic(float *data, int ndata, float *d, int npoles,
					float *future, int nfut);
	float probks(float alam);
	void psdes(unsigned long *lword, unsigned long *irword);
	void  pzextr(int iest, float xest, float *yest, float *yz, float *dy,
					int nv, int nuse);
	void  qcksrt(int n, float *arr);
	float qgaus(float (*func)(float), float a, float b);
	float qromb(float (*func)(float), float a, float b);
	float qromo(float (*func)(float), float a, float b,
					float (*choose)(float (*)(float),float,float,int));
	void  qroot(float *p, int n, float *b, float *c, float eps);
	float qsimp(float (*func)(float), float a, float b);
	float qtrap(float (*func)(float), float a, float b);
	float quad3d(float (*func)(float,float,float), float x1, float x2);
	float ran0(long *idum);
	float ran1(long *idum);
	float ran2(long *idum);
	float ran3(long *idum);
	float ran4(long *idum);
	void  rank(int n, int *indx, int *irank);
	void  ratint(float *xa, float *ya, int n, float x, float *y, float *dy);
	void  realft(float *data, int n, int isign);
	void  red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
					int ic1, int jc1, int jcf, int kc, float ***c, float **s);
	void  rk4(float *y, float *dydx, int n, float x, float h, float *yout,
					void (*derivs)(float,float *,float *));
	void  rkdumb(float *vstart, int nvar, float x1, float x2, int nstep,
					void (*derivs)(float,float *,float *));
	void  rkqc(float *y, float *dydx, int n, float *x, float htry,
					float eps, float *yscal, float *hdid, float *hnext,
					void (*derivs)(float,float *,float *));
	float rofunc(float b);
	float rtbis(float (*func)(float), float x1, float x2, float xacc);
	float rtflsp(float (*func)(float), float x1, float x2, float xacc);
	float rtnewt(void (*funcd)(float,float *,float *), float x1, float x2,
					float xacc);
	float rtsafe(void (*funcd)(float,float *,float *), float x1, float x2,
					float xacc);
	float rtsec(float (*func)(float), float x1, float x2, float xacc);
	void  rzextr(int iest, float xest, float *yest, float *yz, float *dy,
					int nv, int nuse);
	void  scrsho(float (*fx)(float));
	void  shell(int n, float *arr);
	void  shoot(int nvar, float *v, float *delv, int n2, float x1, float x2,
					float eps, float h1, float hmin, float *f, float *dv);
	void  shootf(int nvar, float *v1, float *v2, float *delv1, float *delv2,
					int n1, int n2, float x1, float x2, float xf, float eps,
					float h1, float hmin, float *f, float *dv1, float *dv2);
	void  simp1(float **a, int mm, int *ll, int nll, int iabf, int *kp,
					float *bmax);
	void  simp2(float **a, int n, int *l2, int nl2, int *ip, int kp,
					float *q1);
	void  simp3(float **a, int i1, int k1, int ip, int kp);
	void  simplx(float **a, int m, int n, int m1, int m2, int m3,
					int *icase, int *izrov, int *iposv);
	void  sinft(float *y, int n);
	void  smooft(float *y, int n, float pts);
	void  sncndn(float uu, float emmc, float *sn, float *cn, float *dn);
	void  solvde(int itmax, float conv, float slowc, float *scalv,
					int *indexv, int ne, int nb, int m, float **y, float ***c,
					float **s);
	void  sor(double **a, double **b, double **c, double **d, double **e,
					double **f, double **u, int jmax, double rjac);
	void  sort(int n, float *ra);
	void  sort2(int n, float *ra, float *rb);
	void  sort3(int n, float *ra, float *rb, float *rc);
	void  sparse(float *b, int n, float *x, float *rsq);
	void  spctrm(FILE *fp, float *p, int m, int k, int ovrlap);
	void  spear(float *data1, float *data2, int n, float *d, float *zd,
					float *probd, float *rs, float *probrs);
	void  splie2(float *x1a, float *x2a, float **ya, int m, int n,
					float **y2a);
	void  splin2(float *x1a, float *x2a, float **ya, float **y2a, int m,
					int n, float x1, float x2, float *y);
	void  spline(float *x, float *y, int n, float yp1, float ypn, float *y2);
	void  splint(float *xa, float *ya, float *y2a, int n, float x, float *y);
	void  svbksb(float **u, float *w, float **v, int m, int n, float *b,
					float *x);
	void  svdcmp(float **a, int m, int n, float *w, float **v);
	void  svdfit(float *x, float *y, float *sig, int ndata, float *a,
					int ma, float **u, float **v, float *w, float *chisq,
					void (*funcs)(float,float *,int));
	void  svdvar(float **v, int ma, float *w, float **cvm);
	void  toeplz(float *r, float *x, float *y, int n);
	void  tptest(float *data1, float *data2, int n, float *t, float *prob);
	void  tqli(float *d, float *e, int n, float **z);
	float trapzd(float (*func)(float), float a, float b, int n);
	void  tred2(float **a, int n, float *d, float *e);
	void  tridag(float *a, float *b, float *c, float *r, float *u, int n);
	void  ttest(float *data1, int n1, float *data2, int n2, float *t,
					float *prob);
	void  tutest(float *data1, int n1, float *data2, int n2, float *t,
					float *prob);
	void  twofft(float *data1, float *data2, float *fft1, float *fft2,
					int n);
	void  vander(float *x, float *w, float *q, int n);
	int   zbrac(float (*func)(float), float *x1, float *x2);
	void  zbrak(float (*fx)(float), float x1, float x2, int n, float *xb1,
					float *xb2, int *nb);
	float zbrent(float (*func)(float), float x1, float x2, float tol);
	void  zroots(fcomplex *a, int m, fcomplex *roots, int polish);
