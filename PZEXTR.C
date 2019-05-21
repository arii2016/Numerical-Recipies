#include "nrutil.h"

extern float **d,*x;	/* defined in bsstep */

void pzextr(int iest, float xest, float yest[], float yz[], float dy[],
    int nv, int nuse)
{
	int m1,k1,j;
	float q,f2,f1,delta,*c;

	c=vector(1,nv);
	x[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
		for (j=1;j<=nv;j++)
			d[j][1]=yest[j];
	} else {
		m1=(iest < nuse ? iest : nuse);
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<=m1-1;k1++) {
			delta=1.0/(x[iest-k1]-xest);
			f1=xest*delta;
			f2=x[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=d[j][k1];
				d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) d[j][m1]=dy[j];
	}
	free_vector(c,1,nv);
}
