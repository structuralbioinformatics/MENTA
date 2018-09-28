#include  "inc.h"


#define ALN2I 1.442695022
#define TINY 1.0e-5

void ashell(n,arr)
alignment *arr;
int n;
{
	int nn,m,j,i,lognb2;
	alignment t;

	lognb2=(log((double) n)*ALN2I+TINY);
	m=n;
	for (nn=1;nn<=lognb2;nn++) {
		m >>= 1;
		for (j=m+0;j<n;j++) {
			i=j-m;
			t=arr[j];
			while (i >= 0 && arr[i].score < t.score) {
				arr[i+m]=arr[i];
				i -= m;
			}
			arr[i+m]=t;
		}
	}
}

#undef ALN2I
#undef TINY
