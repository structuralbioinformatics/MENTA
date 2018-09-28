#!/bin/perl
#
 use Getopt::Long;



 &GetOptions(
             "f=s" => \$fasta,
             "m=s" => \$menta,
             "h"   => \$help,
             "o=s" => \$out       );

 if (defined $help){
  printf " 
           f  fasta file with gaps
           m  menta file
           o  output root name
           h  print help\n\n";
  exit();
 }


 open FASTA,"<$fasta";
  while (<FASTA>){
   chop;
   if (/^>/){$name=$_;next;}
   @_=split "",$_;
   for $aa (@_){push @aa,$aa;}
  }
 close FASTA;
 
 open MENTA,"<$menta";
  while (<MENTA>){
   chop;
   if (/^>/){
      push @title,$_;
      $i=0;$j=0;
      next;
     }
   if (/^#/){
      push @title,$_;
      $i++;$j=0;
      next;
     }
   @_=split "",$_;
   for $c (@_){$seq[$i][$j]=$c;$j++;}
  }
 close MENTA;
 open OUT ,">$out";
 for $n (0..$i){
   print OUT "$title[$n]\n";
   $k=0;
   $kk=0;
   for $m (0..$#aa){
       if ($aa[$m] eq "-"){
          print OUT "-";
       }else{
          print OUT "$seq[$n][$k]";
          $k++;
       }
       $kk++;
       if ($kk > 79){print OUT "\n";$kk=0;}
      }
      print OUT "\n";
  }
  close OUT;

