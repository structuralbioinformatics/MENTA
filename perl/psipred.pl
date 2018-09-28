
  $begin=0;
  open (INP,"<$ARGV[0]");
    while (<INP>){
     @data=split;
      
     if($data[0] eq "AA:" && $begin eq 0) {$begin =1;next;}
     if($data[1] eq "PSIPRED" && $begin eq 0) {$begin=1;$psipred=1;next;}
     if ($begin eq 1){
      if($data[0] eq "AA:" ){
       for $i (0..length($data[1])){
         push @aa,substr($data[1],$i,1);
        }
      }
      if($data[0] eq "Pred:"){
       $_=$data[1];
       tr/C/-/;
       for $i (0..length($_)){
         push @ss,substr($_,$i,1);
        }
      }
      if ($psipred eq 1){$psipred++;} 
      if ($psipred eq 2) {
        push @aa,$data[1];
        push @ss,$data[2]; 
      }
     }
    }
  close(INP);

  print ">P1;$ARGV[0]Seq \n";
  print "$ARGV[0]Seq  \n";
  $n=0;
  for $i (0..$#aa){
   print "$aa[$i]";
   $n++;
   if ($n >79 ){print "\n";$n=0;}
  }
  print "*\n";

  print ">P1;$ARGV[0]SS \n";
  print "$ARGV[0]SS  \n";
  $n=0;
  for $i (0..$#ss){
   print "$ss[$i]";
   $n++;
   if ($n >79 ){print "\n";$n=0;}
  }
  print "*\n";


