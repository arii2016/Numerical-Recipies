#include "nrutil.h"

extern float **d,*x;	/* defined in bsstep */

void rzextr(int iest, float xest, float yest[], float yz[], float dy[],
    int nv, int nuse)
{
	int m1,k,j;
	float yy,v,ddy,c,b1,b,*fx;

	fx=vector(1,nuse);
	x[iest]=xest;
	if (iest == 1)
		for (j=1;j<=nv;j++) {
			yz[j]=yest[j];
			d[j][1]=yest[j];
			dy[j]=yest[j];
		}
	else {
		m1=(iest < nuse ? iest : nuse);
		for (k=1;k<=m1-1;k++)
			fx[k+1]=x[iest-k]/xest;
		for (j=1;j<=nv;j++) {
			v=d[j][1];
			d[j][1]=yy=c=yest[j];
			for (k=2;k<=m1;k++) {
				b1=fx[k]*v;
				b=b1-c;
				if (b) {
					b=(c-v)/b;
					ddy=c*b;
					c=b1*b;
				} else
					ddy=v;
				if (k != m1) v=d[j][k];
				d[j][k]=ddy;
				yy += ddy;
			}
			dy[j]=ddy;
			yz[j]=yy;
		}
	}
	free_vector(fx,1,nuse);
}
