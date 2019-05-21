extern float **a;

void asub(float x[], float v[], int n)
{
	int i,j;
	for (i=1;i<=n;i++) {
		v[i]=0.0;
		for (j=1;j<=n;j++) v[i] += a[i][j]*x[j];
	}
}
