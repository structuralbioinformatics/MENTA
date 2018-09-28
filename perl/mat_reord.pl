#!/usr/local/bin/perl

open(MAT,"<$ARGV[0]") || die "Cannot open $ARGV[0].$!\n";
$newfile="$ARGV[0]".".reord";
open (NEWMAT,">$newfile");

print "Want to write symbols in reordered matrix ? ";
$symwrt=<STDIN>; chop $symwrt; 
print "New symbols:\n";
$symb=<STDIN>; chop $symb;

print NEWMAT "\# Reordered matrix from $ARGV[0]\n";

while (<MAT>)  {
	if ($_=~ /^\#/)  { print NEWMAT "$_"; next; }
	if ($_!~ /\d/)  { &fILE_order($_);next; }
	$_=~ s/-/ -/g;
	$_=~ s/[A-Z]//g;
	push(@Mat,[split " ",$_]);
	}
close MAT;

  
for $i(0..$#ordsym)  {
	if ($ordsym[$i]=~ / /)  { next; }
	if ($symwrt eq "y")  { print NEWMAT "$ordsym[$i] "; }
	for $j(0..$#ordsym)  {
		if ($ordsym[$j]=~ /\s/)  { next; }
		$m=$hash{$ordsym[$i]};
		$n=$hash{$ordsym[$j]};
                if ($m >= 0 && $n >= 0 ){
                 if ($Mat[$n][$m] eq 0.0 && $Mat[$m][$n] ne 0.0){ $Mat[$n][$m]=$Mat[$m][$n];}
                 if ($Mat[$m][$n] eq 0.0 && $Mat[$n][$m] ne 0.0){ $Mat[$m][$n]=$Mat[$n][$m];}
                 if ($Mat[$m][$n] ne 0.0 && $Mat[$n][$m] ne 0.0){ $x=($Mat[$m][$n]+$Mat[$n][$m])/2.0;
                                                                  $Mat[$m][$n]=$x;
                                                                  $Mat[$n][$m]=$x;}
		 printf NEWMAT ("%10.2f",$Mat[$n][$m]);
                }elsif ($ordsym[$i] eq "-" || $ordsym[$j] eq "-"){
                 printf NEWMAT ("%10.2f",50.0);
                }else{
                 printf NEWMAT ("%10.2f",0.0);
                }
		}
	print NEWMAT "\n";
	}
if ($symwrt eq "y")  {
	foreach $ordsym(@ordsym)  {  print NEWMAT "    $ordsym"; }
	}
close NEWMAT;


### Subroutines from now ##############################
sub fILE_order
{
$_=~ s/ //g;
@symbols=split("",$_);
@ordsym=split("",$symb);
print "\nReordered symbols:\n";
print "@ordsym\n";

foreach $ordsym(@ordsym)  {
        $hash{$ordsym}=-1;
	for $n(0..$#symbols)  {
		if ($symbols[$n] eq $ordsym)  { 
			$hash{$ordsym}=$n; 
			#print "$ordsym $symbols $hash{$ordsym} $n\n";
			next; }
		}
          #print "$ordsym $hash{$ordsym} \n";
	}
}


### By Wisl, Jul 2001. 

