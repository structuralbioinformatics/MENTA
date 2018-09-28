#!/usr/bin/perl
#

 use Getopt::Long;



 &GetOptions("n=s" => \$groups,
             "d=s" => \$interval,
             "h"   => \$help,
             "o=s" => \$output       );

 if (defined $help){
  printf " 
           n  Number of methods to Use in MENTA
           d  Interval of weights for each method, from 0 to 100.
             (ie. 50 means there is a combination with a weight
             for a method double than for the rest)  
           o  output root name
           h  print help\n\n";
  exit();
 }



 undef @set;
 $j=0;
 $total=1;
 undef %input;

 for $i (0..$total){$suffix[$i]="";}
 $set[0]=0;
 while ($j<100){ $j+=$interval;push @set,$j;}
 for $i (1..$groups){
      $k=0;
      for ($i=0;$i<$total;$i++){
      for ($j=0;$j<$#set+1;$j++){
        $suffixn[$k]="$suffix[$i]";
        $k++;
      }
      }
      $k=0;
      for ($i=0;$i<$total;$i++){
      for ($j=0;$j<$#set+1;$j++){
        $suffix[$k]="$suffixn[$k]"."-"."$set[$j]";
        $k++;
      }
      }
     $total=$k;
 }

 for $word (@suffix){
  undef @combo;
  @combo=split /-/,$word;
  $sum=0;
  for $x (1..$#combo){ $sum+=$combo[$x];}
  if ($sum>0){
   undef @combi;
   for $x (1..$#combo){ $y=int(100*($combo[$x])/$sum);push @combi,$y;}
   $string="";
   for $y (@combi){ $string.="-"."$y"; }
   $input{$string}=1;
  }
 }
 undef @suffix;
 $i=0;
 foreach $key (keys %input){
  $i++;
  @combo=split /-/,$key;
  $out="$output"."_"."$i".".out";
  printf " INPUT $i ";
  for  $x (1..$#combo){printf "  $combo[$x] "; }
  printf "\n";
  open (OUT, ">$out");
    for  $x (1..$#combo){printf OUT " $combo[$x] \n";}
  close OUT;
 }

