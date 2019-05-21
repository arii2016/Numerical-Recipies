#include <math.h>
#include "complex.h"

#define NPMAX 100
#define ZERO Complex(0.0,0.0)
#define ONE Complex(1.0,0.0)

void fixrts(float d[], int npoles)
{
	void zroots(fcomplex a[], int m, fcomplex roots[], int polish);
	int i,j,polish;
	fcomplex a[NPMAX],roots[NPMAX];

	a[npoles]=ONE;
	for (j=npoles-1;j>=0;j--)
		a[j]=Complex(-d[npoles-j],0.0);
	polish=1;
	zroots(a,npoles,roots,polish);
	for (j=1;j<=npoles;j++)
		if (Cabs(roots[j]) > 1.0)
			roots[j]=Cdiv(ONE,Conjg(roots[j]));
	a[0]=Csub(ZERO,roots[1]);
	a[1]=ONE;
	for (j=2;j<=npoles;j++) {
		a[j]=ONE;
		for (i=j;i>=2;i--)
			a[i-1]=Csub(a[i-2],Cmul(roots[j],a[i-1]));
		a[0]=Csub(ZERO,Cmul(roots[j],a[0]));
	}
	for (j=0;j<=npoles-1;j++)
		d[npoles-j] = -a[j].r;
}
