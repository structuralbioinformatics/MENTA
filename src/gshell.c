#include  "inc.h"


#define ALN2I 1.442695022
#define TINY 1.0e-5

void gShell(ii,n,g,p)
int ii;
int *n,**g;
profile **p;
{
	int nn,m,j,i,lognb2;
	int t;



	lognb2=(log((double) n[ii])*ALN2I+TINY);
	m=n[ii];
	for (nn=1;nn<=lognb2;nn++) {
		m >>= 1;
		for (j=m+0;j<n[ii];j++) {
			i=j-m;
			t=g[ii][j];
			//while (    i >= 0 && (p[g[ii][i]][0].score<p[t][0].score) ) {
			while (    i >= 0 && (p[g[ii][i]][0].ident < p[t][0].ident || (p[g[ii][i]][0].ident==p[t][0].ident && p[g[ii][i]][0].homol<p[t][0].homol) || (p[g[ii][i]][0].ident==p[t][0].ident && p[g[ii][i]][0].homol==p[t][0].homol && p[g[ii][i]][0].score<p[t][0].score)) ) {
			//ºwhile (i >= 0 && p[g[ii][i]][0].ident < p[t][0].ident ) {
				g[ii][i+m]=g[ii][i];
				i -= m;
			}
			g[ii][i+m]=t;
		}
	}
}

#undef ALN2I
#undef TINY
