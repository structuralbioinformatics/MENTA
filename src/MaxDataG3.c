#include  "inc.h"

max3  MaxDataG3(data,limit,tolerance)
float   data[6];
float   limit;
float   *tolerance;
{
 float  max;
 max3   maxim;
 float  factor,factor0,factor1;
 int    d0,d1,d2,imax;

 factor   = tolerance[0];
 factor0  = pow( 10,factor);
 factor1  = pow( 10,tolerance[1]);

 
 max=data[0];
 if (data[1]>max) {max=data[1];}
 if (data[2]>max) {max=data[2];}
 if (max>limit){max=limit;} 
 maxim.value=max;


 if (max>1.0e-20){
  imax=(int)max;
  while(imax<10){
   factor += 1.0 ;
   factor0 = pow( 10,factor);
   imax=(int)(max*factor0);
  }
  factor += tolerance[0] ;
  factor0 = pow( 10,factor);
  imax=(int)(max*factor0);
 }


 d0   =  (int)( data[0]*factor0 );
 d1   =  (int)( data[1]*factor0 );
 d2   =  (int)( data[2]*factor0 );

   imax = d0;

   if ( d1 > imax) {imax=d1;}
   if ( d2 > imax) {imax=d2;}
   if ( imax-d0 < (int)(factor1) ) {maxim.first=1;}else{maxim.first=0;}
   if ( imax-d1 < (int)(factor1) ) {maxim.second=1;}else{maxim.second=0;}
   if ( imax-d2 < (int)(factor1) ) {maxim.third=1;}else{maxim.third=0;}

 max=data[3];
 if (data[3]>max) {max=data[3];}
 if (data[4]>max) {max=data[4];}
 if (data[5]>max) {max=data[5];}
 maxim.score=max;

 if (imax==0){
    maxim.first=1;
    maxim.second=0;
    maxim.third=0;
 }else{
    if (maxim.first==0 && maxim.second==0 && maxim.third==0){
     if (imax==d0){ maxim.first=1;}
     if (imax==d1){ maxim.second=1;}
     if (imax==d2){ maxim.third=1;}
    }
 }

  if (tolerance[4] >= 1.0) {

    if (maxim.first==1 && maxim.second==1 && maxim.third==1){
     if (imax==d0){ maxim.first=1;}
     if (imax==d1){ maxim.second=1;}
     if (imax==d2){ maxim.third=1;}
    }
    if (maxim.first==1 && maxim.second==1 && maxim.third==1){
         maxim.first=1;
	 maxim.second=0;
	 maxim.third=0;
    }

  if (tolerance[4] >= 2.0 ){

    if (maxim.second==1 && maxim.third==1){
        if (imax==d1){ maxim.second=1;}
        if (imax==d2){ maxim.third=1;}
    }
    if (maxim.second==1 && maxim.third==1){
	 maxim.first=0;
	 maxim.second=1;
	 maxim.third=0;
    }

  if (tolerance[4] >= 3.0) {

     if (imax==d0){ maxim.first=1;}
     if (imax==d1){ maxim.second=1;}
     if (imax==d2){ maxim.third=1;}

    if (maxim.first==1 && maxim.third==1){
	 maxim.first=1;
	 maxim.second=0;
	 maxim.third=0;
    }
    if (maxim.first==1 && maxim.second==1){
	 maxim.first=1;
	 maxim.second=0;
	 maxim.third=0;
    }

  }}}


 return maxim;

}


