#!/bin/perl
#
 use Getopt::Long;



 &GetOptions("i=s" => \$input ,
             "f=s" => \$name,
             "m=s" => \$matrix_m,
             "s=s" => \$matrix_s,
             "b=s" => \$matrix_b,
             "w=s" => \$weight,
             "h"   => \$help,
             "o=s" => \$out       );

 if (defined $help){
  printf " 
           i  INPUT FASTA file
           f  Working name to tag the files
           m  Mizuguchi matrix
           b  Buriedness matrix
           s  Secondary Structure matrix
           w  weight/s
           o  output root name
           h  print help\n\n";
  exit();
 }

 $sol={};
 $sol={ A => 115.0, R => 263.0, N => 184.0, D => 170.0, C => 149.0, 
        Q => 208.0, E => 207.0,
        G => 86.0, H => 206.0, I => 187.0, L => 192.0, K => 222.0, 
        M => 210.0, F => 230.0,
        P => 140.0 , S => 140.0, T => 164.0, W => 269.0, Y => 257.0, V => 161.0};
 $bur={};
 $bur={ 9 => 'A', 8 => 'B', 7 => 'C', 6 => 'D', 5 => 'E', 
        4 => 'F', 3 => 'G', 2 => 'H', 1 => 'I', 0 => 'W'};

 if (defined $ENV{SABLE}){$sable="$ENV{SABLE}"."/run.sable";}else{$sable="/home/boliva/PROGRAMS/Sequence/SABLE/sable_distr";}
 if (defined $ENV{MENTA}){$home="$ENV{MENTA}";}else{$home="/home/boliva/PROJECTS/MENTA_M/MMEntA/MEntA";}
 $sable2menta="$home"."/perl/SABLE2Menta.pl";

 undef %fasta;
 undef @prot;
 open(INP,"<$input");
     undef @seq;
     $seqid=" ";
     $n=0;
     while(<INP>){
      if (/^>(\S+)/){
	      $fasta{$seqid}=join "",@seq;
	      $seqid=$1;
              push @prot,$seqid;
              printf ".";
	      undef @seq;
              $n++;
	      next;}
      chop;
      push @seq,$_;
     }
  $fasta{$seqid}=join "",@seq;
 close INP;
 for ($i=0;$i<$n;$i++){
     $id=$prot[$i];
     $protein="$name"."."."$i".".fa";
     $outsable="$name"."."."$i".".faSABLE.out";
     $outmint="$name"."."."$id";
     $outmenta="$name"."."."$id".".menta";
     printf "\nGenerate Fasta file $protein ID $id\n";
     open (OUT,">$protein");
     printf OUT ">$id\n";
     $j=0;
     @residue=split ("",$fasta{$id});
     for $res (@residue){
              if ($j>79){printf OUT "\n";$j=0}
              printf OUT "$res";
              $j++;
          }
          if ($#residue>0) {printf OUT "\n";}
     close OUT;
     printf "Run Sable on $protein\n";
     system("\\rm -r OUT_SABLE");
     if (-e $outsable){printf "Use $outsable\n";}else{system("$sable $protein");}
     printf "Run Sable2Menta -i  $outsable -o $outmint \n";
     $matrix=" ";
     if (defined $matrix_m) {$matrix.=" -m $matrix_m";}
     if (defined $matrix_s) {$matrix.=" -b $matrix_b";}
     if (defined $matrix_b) {$matrix.=" -s $matrix_s";}
     if ($weight<=0){ $weight=0;}
     printf "perl $sable2menta -i  $outsable -f $id $matrix -w $weight  -o  $outmint";
     if (-e $outmenta){printf "Use $outmenta\n";}else{system("perl $sable2menta -i  $outsable -f $id $matrix -w $weight  -o  $outmint");}
     if ($i==0){system("cat $outmenta > $out");}
     else      {system("cat $outmenta >> $out");}
 }
  

