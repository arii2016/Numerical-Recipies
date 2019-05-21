#include <math.h>
#include "nrutil.h"
#define IMAX 11
#define NUSE 7
#define SHRINK 0.95
#define GROW 1.2

float **d=0,*x=0;	/* defining declaration */

void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []))
{
	void mmid(float y[], float dydx[], int nvar, float xs, float htot,
		int nstep, float yout[], void (*derivs)(float, float[], float[]));
	void rzextr(int iest, float xest, float yest[], float yz[], float dy[],
		int nv, int nuse);
	int i,j;
	float xsav,xest,h,errmax,temp;
	float *ysav,*dysav,*yseq,*yerr;
	static int nseq[IMAX+1]={0,2,4,6,8,12,16,24,32,48,64,96};

	ysav=vector(1,nv);
	dysav=vector(1,nv);
	yseq=vector(1,nv);
	yerr=vector(1,nv);
	x=vector(1,IMAX);
	d=matrix(1,nv,1,NUSE);
	h=htry;
	xsav=(*xx);
	for (i=1;i<=nv;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	for (;;) {
		for (i=1;i<=IMAX;i++) {
			mmid(ysav,dysav,nv,xsav,h,nseq[i],yseq,derivs);
			xest=(temp=h/nseq[i],temp*temp);
			rzextr(i,xest,yseq,y,yerr,nv,NUSE);
			if (i > 3) {
				errmax=0.0;
				for (j=1;j<=nv;j++)
					if (errmax < fabs(yerr[j]/yscal[j]))
						errmax=fabs(yerr[j]/yscal[j]);
				errmax /= eps;
				if (errmax < 1.0) {
					*xx += h;
					*hdid=h;
					*hnext = i==NUSE? h*SHRINK : i==NUSE-1?
						h*GROW : (h*nseq[NUSE-1])/nseq[i];
					free_matrix(d,1,nv,1,NUSE);
					free_vector(x,1,IMAX);
					free_vector(yerr,1,nv);
					free_vector(yseq,1,nv);
					free_vector(dysav,1,nv);
					free_vector(ysav,1,nv);
					return;
				}
			}
		}
		h *= 0.25;
		for (i=1;i<=(IMAX-NUSE)/2;i++) h /= 2.0;
		if ((float)(*xx+h) == (*xx)) nrerror("Step size underflow in bsstep");
	}
}

#undef IMAX
#undef NUSE
#undef SHRINK
#undef GROW
