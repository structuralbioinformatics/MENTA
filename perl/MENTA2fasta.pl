#!/bin/perl
#
 use Getopt::Long;



 &GetOptions("i=s" => \$input ,
             "m=s" => \$matrix_m,
             "s=s" => \$matrix_s,
             "b=s" => \$matrix_b,
             "w=s" => \$weight,
             "h"   => \$help,
             "o=s" => \$out       );

 if (defined $help){
  printf " 
           i  INPUT MENTA file
           m  Mizuguchi matrix (1st)
           s  Secondary Structure matrix (2on)
           b  Buriedness matrix (3rd)
           w  weight/s
           o  output root name
           h  print help\n\n";
  exit();
 }


 open(INP,"<$input");
 open(OUT,">$out");
 $found=0;
 $n=0;
 undef %extra;
 undef %start;
 undef %profile;
 while(<INP>){
   if (/^# PROFILE \[(\d+)\]/){
      printf "Find profile $1\n";
      if ($id[$#id] eq $1) {
       if ($n==1) {$extra{m}="# $weight  $matrix_m";next;}
       if ($n==2) {$extra{s}="# $weight  $matrix_s";next;}
       if ($n==3) {$extra{b}="# $weight  $matrix_b";next;}
      }else{ 
       push @id,$1;
       $n=0;
       $found=1;
       foreach $p (keys %start){
       printf OUT "$start{$p}{seq}\n";
       $j=0;
       @residue=split ("",$profile{$p}{seq});
       for $res (@residue){
              if ($j>79){printf OUT "\n";$j=0}
              printf OUT "$res";
              $j++;
          }
          if ($#residue>0) {printf OUT "\n";}
       foreach $type (%extra){
         undef @residue;
         if (defined $extra{$type}){
          printf OUT "$extra{$type}\n";
          $j=0;
          @residue=split ("",$profile{$p}{m});
          for $res (@residue){
              if ($j>79){printf OUT "\n";$j=0}
              printf OUT "$res";
              $j++;
           }
           if ($#residue>0) {printf OUT "\n";}
          }
         }
       }
       undef %extra;
       undef %start;
       undef %profile;
      }
      next;
   }
   if (/^\*\*\*\*\*\*\*\*/){$n++;$found=0;}
   if ($found==1 && $n==0 ){
     @line=split;
     if ($line[0] ne "" && $line[0] ne "CONSENSUS_____"){
        $start{$line[0]}{seq} = ">M$id[$#id]"." "."$line[0]";
        $profile{$line[0]}{seq}.=$line[3];
     }
   }
   if ($found==0 && $n==1){
     @line=split;
     if ($line[0] ne "" && $line[0] ne "CONSENSUS_____"){
        $profile{$line[0]}{m}.=$line[3];
     }
   }
   if ($found==0 && $n==2){
     @line=split;
     if ($line[0] ne "" && $line[0] ne "CONSENSUS_____"){
        $profile{$line[0]}{s}.=$line[3];
     }
   }
   if ($found==0 && $n==3){
     @line=split;
     if ($line[0] ne "" && $line[0] ne "CONSENSUS_____"){
        $profile{$line[0]}{b}.=$line[3];
     }
   }
 }

 foreach $p (keys %start){
    printf OUT "$start{$p}{seq}\n";
    $j=0;
    @residue=split ("",$profile{$p}{seq});
    for $res (@residue){
              if ($j>79){printf OUT "\n";$j=0}
              printf OUT "$res";
              $j++;
       }
       if ($#residue>0) {printf OUT "\n";}
    foreach $type (%extra){
      undef @residue;
      if (defined $extra{$type}){
       printf OUT "$extra{$type}\n";
       $j=0;
       @residue=split ("",$profile{$p}{m});
       for $res (@residue){
           if ($j>79){printf OUT "\n";$j=0}
           printf OUT "$res";
           $j++;
        }
        if ($#residue>0) {printf OUT "\n";}
       }
      }
    }
    undef %extra;
    undef %start;
    undef %profile;


