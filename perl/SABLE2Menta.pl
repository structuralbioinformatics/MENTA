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
           i  INPUT SABLE file
           f  Sequence name
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


  $read_ss=0;
  $read_sa=0;
  undef @aa,@ss,@bb;
  open (INP,"$input");
    while (<INP>){
     @data=split;
     if ($data[0] eq "SECTION_SS"){$read_ss=1;} 
     if ($data[0] eq "SECTION_SA"){$read_sa=1;} 
     if ($data[0] eq "END_SECTION"){$read_ss=$read_sa=0; }
     if (substr($data[0],0,1) eq ">" && ($read_ss eq 1 || $read_sa eq 1)){$begin=1;next;}
     if ($begin eq 2 ){
       for $i (0..length($data[0])){
         if ($read_ss eq 1 && substr($data[0],$i,1) ne ""){ push @ss,substr($data[0],$i,1);}
         if ($read_sa eq 1 && substr($data[0],$i,1) ne ""){ push @bb,substr($data[0],$i,1);}
        }
       $begin=0;
       next;
      }
     if ($begin eq 1 ){
       for $i (0..length($data[0])){
         if ($read_ss eq 1 && substr($data[0],$i,1) ne ""){push @aa,substr($data[0],$i,1);}
        }
       $begin=2;
       next;
      }
    }
  close(INP);

  $data={};
  $data={ AA  => [@aa],
          SS  => [@ss],
          ACC => [@bb],};


  (@mizu) = &Mizuguchi($data);

  $outpir="$out".".pir";
  $outmenta="$out".".menta";

  open OUTP, ">$outpir";

  print OUTP ">P1;$name\n";
  print OUTP "$name Seq \n";
  $n=0;
  for $i (0..$#aa){
   print OUTP "$aa[$i]";
   $n++;
   if ($n >79 ){print OUTP "\n";$n=0;}
  }
  print OUTP "*\n";

  print OUTP ">P1;$name SS\n";
  print OUTP "$name SS \n";
  $n=0;
  for $i (0..$#ss){
   print OUTP "$ss[$i]";
   $n++;
   if ($n >79 ){print OUTP "\n";$n=0;}
  }
  print OUTP "*\n";
  
  print OUTP ">P1;$name BUR\n";
  print OUTP "$name BUR \n";
  $n=0;
  for $i (0..$#bb) {
   print OUTP "$bur->{$bb[$i]}";
   $n++;
   if ($n >79 ){print OUTP "\n";$n=0;}
  }
  print OUTP "*\n";

  print OUTP ">P1;$name ACC\n";
  print OUTP "$name ACC \n";
  $n=0;
  for $i (0..$#bb) {
   print OUTP "$bb[$i]";
   $n++;
   if ($n >79 ){print OUTP "\n";$n=0;}
  }
  print OUTP "*\n";
 
 
  print OUTP ">P1;$name SS_ACC\n";
  print OUTP "$name SS_ACC \n";
  $n=0;
  for $i (0..$#mizu) {
   print OUTP "$mizu[$i]";
   $n++;
   if ($n >79 ){print OUTP "\n";$n=0;}
  }
  print OUTP "*\n";
  
  close OUTP; 


  open OUTM, ">$outmenta";
  print OUTM ">$name\n";
  $n=0;
  for $i (0..$#aa){
   print OUTM "$aa[$i]";
   $n++;
   if ($n >79 ){print OUTM "\n";$n=0;}
  }
  print OUTM "\n";
  if (defined $matrix_s){
   printf OUTM "#%8.5f  %s\n",$weight,$matrix_s;
   $n=0;
   for $i (0..$#ss) {
    print OUTM "$ss[$i]";
    $n++;
    if ($n >79 ){print OUTM "\n";$n=0;}
   }
   print OUTM "\n";
  }
  if (defined $matrix_b){
   printf OUTM "#%8.5f  %s\n",$weight,$matrix_b;
   $n=0;
   for $i (0..$#bb) {
    print OUTM "$bur->{$bb[$i]}";
    $n++;
    if ($n >79 ){print OUTM "\n";$n=0;}
   }
   print OUTM "\n";
  }
  if (defined $matrix_m){
   printf OUTM "#%8.5f  %s\n",$weight,$matrix_m;
   $n=0;
   for $i (0..$#mizu) {
    print OUTM "$mizu[$i]";
    $n++;
    if ($n >79 ){print OUTM "\n";$n=0;}
   }
   print OUTM "\n";
  }

  close OUTM;


sub Mizuguchi
{
 my $data=$_[0];
 my ($i,$j,$h310,$strand,$helix,$alfa,$beta,@mizu,@aa,@bb,@ss);


  (@ss)=@{$data->{SS}};
  (@bb)=@{$data->{ACC}};
  (@aa)=@{$data->{AA}};

  $h310=0;
  $strand=0;
  $helix=0;
  $beta=0;
  $alfa=0;
  $end=0;
  undef @mizu;
  for $i (0..$#ss){
    #printf "Test $i $ss[$i] Alfa $alfa Beta $beta Strand $strand Helix $helix \n";
    #for ($ii=0;$ii<$i;$ii++){printf "$ss[$ii]";}
    #printf "\n";
    #for ($ii=0;$ii<=$#mizu;$ii++){printf "$mizu[$ii]";}
    #printf "\n";
    if ($ss[$i] eq "G") {$h310=1;$strand=0;$helix=0;}
    if ($ss[$i] eq "E") {$h310=0;$strand=1;$helix=0;}
    if ($ss[$i] eq "H") {$h310=0;$strand=0;$helix=1;}
    if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {
        $h310=0;$strand=0;$helix=0;
        if ($alfa==0 && $beta == 0) {push @mizu,"C";if ($i eq $#ss){$end=1;}}
        }
    if ($h310 > 0) {push @mizu,"G";}
    if ($strand > 0 ) { $beta++; }
    if ($helix  > 0 ) { $alfa++; }
    if ($beta == 1) {$b0=$i;$acc_beta=$bb[$i];}
    if ($alfa == 1) {$a0=$i;if ($bb[$i]>5) {$acc_alfa=1;}else{$acc_alfa=0;} }
    if ($strand > 0 && $beta > 1){ $acc_beta+=$bb[$i];   }    
    if ($helix > 0 && $alfa > 1){ if ($bb[$i]>5) {$acc_alfa+=1;} }
    if ($strand == 0 && $beta > 0){
        $size=$i-$b0;
        if ($size < 1) {
         $beta=0;
         push @mizu,"C";
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }else{
         $acc_beta=$acc_beta/$size;
         $short=$medium=$large=0;
         if ($size <= 3){$short=1;}
         if ($size <= 8 && $size > 3){$medium=1;}
         if ($size > 8) {$large=1;}
         if ($short == 1){if ($acc_beta <=2){$type="D";}else{$type="S";}}
         if ($medium== 1){if ($acc_beta <=2){$type="B";}else{$type="E";}}
         if ($large == 1){if ($acc_beta <=2){$type="M";}else{$type="N";}}
         for $j ($b0..($i-1)){push @mizu,$type; }
         $beta=0;
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }
       }
    if ($helix == 0 && $alfa > 0){
        $size=$i-$a0;
        if ($size < 1) {
         $alfa=0;
         push @mizu,"C";
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }else{
         $short=$large=0;
         if ($size <= 7){$short=1;}
         if ($size > 7) {$large=1;}
         if ($short == 1){if ($acc_alfa < 3){$type="A";}else{$type="F";}}
         if ($large == 1){if ($acc_alfa/2 < 3){$type="H";}else{$type="I";}}
         for $j ($a0..($i-1)){push @mizu,$type;}
         $alfa=0;
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }
       }
  }
    $i=$#ss;
    if ($h310 > 0) {push @mizu,"G";}
    if ($strand > 0 ) { $beta++; }
    if ($helix  > 0 ) { $alfa++; }
    if ($beta == 1) {$b0=$i;$acc_beta=$bb[$i];}
    if ($alfa == 1) {$a0=$i;if ($bb[$i]>5) {$acc_alfa=1;}else{$acc_alfa=0;} }
    if ($strand > 0 && $beta > 1){ $acc_beta+=$bb[$i];   }    
    if ($helix > 0 && $alfa > 1){ if ($bb[$i]>5) {$acc_alfa+=1;} }
    if ($beta > 0){
        $size=$i-$b0;
        if ($size < 1) {
         push @mizu,"C";
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }else{
         $acc_beta=$acc_beta/$size;
         $short=$medium=$large=0;
         if ($size <= 3){$short=1;}
         if ($size <= 8 && $size > 3){$medium=1;}
         if ($size > 8) {$large=1;}
         if ($short == 1){if ($acc_beta <=2){$type="D";}else{$type="S";}}
         if ($medium== 1){if ($acc_beta <=2){$type="B";}else{$type="E";}}
         if ($large == 1){if ($acc_beta <=2){$type="M";}else{$type="N";}}
         for $j ($b0..($i)){push @mizu,$type;}
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}else{push @mizu,"$type";}
        }
       }
    if ($alfa > 0){
        $size=$i-$a0;
        if ($size < 1) {
         push @mizu,"C";
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}
        }else{
         $short=$large=0;
         if ($size <= 7){$short=1;}
         if ($size > 7) {$large=1;}
         if ($short == 1){if ($acc_alfa < 3){$type="A";}else{$type="F";}}
         if ($large == 1){if ($acc_alfa/2 < 3){$type="H";}else{$type="I";}}
         for $j ($a0..($i)){push @mizu,$type;}
         if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {push @mizu,"C";}else{push @mizu,"$type";}
        }
       }
     if ($alfa==0 && $beta == 0 && $end == 0) {push @mizu,"C";}
   return (@mizu);
} 

