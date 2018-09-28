#include  "inc.h"

int ConsensusID(x,s)
int   x,s;
{
 int c;

 c=0;

 if      (s>=100)        { if (x>50)c=1; }
 else if (s<100 && s>=50){ if (x>55)c=1; } 
 else if (s<50 && s>=30) { if (x>60)c=1; } 
 else if (s<30 && s>=20) { if (x>65)c=1; } 
 else if (s<20 && s>=10) { if (x>70)c=1; } 
 else if (s<10 && s>=5)  { if (x>75)c=1; } 
 else if (s<5  && s>=1)  { if (x>80)c=1; } 

 return c;
}

