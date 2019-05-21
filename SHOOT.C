#include "nrutil.h"

void shoot(int nvar, float v[], float delv[], int n2, float x1, float x2,
    float eps, float h1, float hmin, float f[], float dv[])
{
	int nok,nbad,iv,i,*indx;
	float sav,det,*y,**dfdv;
	void odeint(float ystart[], int nvar, float x1, float x2,
		float eps, float h1, float hmin, int *nok, int *nbad,
		void (*derivs)(float, float [], float []),
		void (*rkqs)(float [], float [], int, float *, float, float,
		float [], float *, float *, void (*)(float, float [], float [])));
	void ludcmp(float **a, int n, int *indx, float *d);
	void lubksb(float **a, int n, int *indx, float b[]);
	void derivs(float x, float y[], float dydx[]);
	void rkqc(float y[], float dydx[], int n, float *x,
		float htry, float eps, float yscal[], float *hdid, float *hnext,
		void (*derivs)(float, float [], float []));
	void load(float x1, float v[], float y[]);
	void score(float xf, float y[], float f[]);

	y=vector(1,nvar);
	indx=ivector(1,nvar);
	dfdv=matrix(1,n2,1,n2);
	load(x1,v,y);
	odeint(y,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqc);
	score(x2,y,f);
	for (iv=1;iv<=n2;iv++) {
		sav=v[iv];
		v[iv] += delv[iv];
		load(x1,v,y);
		odeint(y,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqc);
		score(x2,y,dv);
		for (i=1;i<=n2;i++)
			dfdv[i][iv]=(dv[i]-f[i])/delv[iv];
		v[iv]=sav;
	}
	for (iv=1;iv<=n2;iv++) dv[iv] = -f[iv];
	ludcmp(dfdv,n2,indx,&det);
	lubksb(dfdv,n2,indx,dv);
	for (iv=1;iv<=n2;iv++) v[iv] += dv[iv];
	free_matrix(dfdv,1,n2,1,n2);
	free_ivector(indx,1,nvar);
	free_vector(y,1,nvar);
}
