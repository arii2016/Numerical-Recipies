void mdian1(float x[], int n, float *xmed)
{
	int n2,n2p;
	void hpsort(unsigned long n, float ra[]);

	hpsort(n,x);
	n2p=(n2=n/2)+1;
	*xmed=(n % 2 ? x[n2p] : 0.5*(x[n2]+x[n2p]));
}
