#include "nrutil.h"

void lfit(float x[], float y[], float sig[], int ndata, float a[],
	int ma, int lista[], int mfit, float **covar, float *chisq,
	void (*funcs)(float,float *,int))
{
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	int k,kk,j,ihit,i;
	float ym,wt,sum,sig2i,**beta,*afunc;

	beta=matrix(1,ma,1,1);
	afunc=vector(1,ma);
	kk=mfit+1;
	for (j=1;j<=ma;j++) {
		ihit=0;
		for (k=1;k<=mfit;k++)
			if (lista[k] == j) ihit++;
		if (ihit == 0)
			lista[kk++]=j;
		else if (ihit > 1) nrerror("Bad lista permutation in LFIT-1");
	}
	if (kk != (ma+1)) nrerror("Bad lista permutation in LFIT-2");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma)
			for (j=(mfit+1);j<=ma;j++)
				ym -= a[lista[j]]*afunc[lista[j]];
		sig2i=1.0/SQR(sig[i]);
		for (j=1;j<=mfit;j++) {
			wt=afunc[lista[j]]*sig2i;
			for (k=1;k<=j;k++)
				covar[j][k] += wt*afunc[lista[k]];
			beta[j][1] += ym*wt;
		}
	}
	if (mfit > 1)
		for (j=2;j<=mfit;j++)
			for (k=1;k<=j-1;k++)
				covar[k][j]=covar[j][k];
	gaussj(covar,mfit,beta,1);
	for (j=1;j<=mfit;j++) a[lista[j]]=beta[j][1];
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,lista,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
}
