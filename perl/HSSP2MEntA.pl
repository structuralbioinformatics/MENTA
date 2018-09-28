#!/usr/bin/perl -w

 use Getopt::Long;



 &GetOptions("i=s" => \$hssp,
             "m=s" => \$matrix_m,
             "s=s" => \$matrix_s,
             "b=s" => \$matrix_b,
             "w=s" => \$weight,
             "h"   => \$help,
             "o=s" => \$out       );

 if (defined $help){
  printf " 
           i  INPUT HSSP file
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


my @alineament = read_HSSP($hssp);

my $i;

for($i=0; $i<scalar(@alineament); $i++){
        $output="$out"."_"."$i".".menta";
	write_align_HSSP(\@alineament, $output, $i);
}

#########   FUNCIO PER LLEGIR UN FITXER EN FORMAT HSSP #############
## Retorna una matriu amb l'alineament
## Par�etres d'entrada: el nom del fitxer HSSP d'origen
sub read_HSSP {

	my $fitxer_align = $_[0];		#nom del fitxer de l'alineament d'origen (primer parametre)
	my $lectura_linia=" ";			#variable per emmagatzemar les lectures de les linies

	#variables per a bucles
	my $i;
	my $j;

	my @alineament;     #matriu on s'emmagatzema tot l'alineament final
	my $posicio=0;      #variable que indicara quina posicio de l'alineament s'est�treballant
	my @linia_camps;    #vector on s'emmagatzema temporalment tota una linia de lectura amb els seus camps separats
	my @secuencia;      #vector on s'emmagatzema temporalment la sequencia d'una linia, separada per lletres
	my @interval_prots; #vector que tindr�dos posicions, on s'indica el m�im i m�im de l'identificador de proteina que es llegir�a cada bloc
        $step=0;
	open(ALINEAMIENTO, $fitxer_align) || die "no es pot obrir el fitxer: $!";

	while( $lectura_linia = <ALINEAMIENTO> ){

	    ## LLEGEIX EL NOM DE LES PROTEINES I ELS POSA A LA PRIMERA FILA DE LA MATRIU DE L'ALINEAMENT
	    if ($lectura_linia =~ m/PROTEINS\s:/ ){
	       $lectura_linia = <ALINEAMIENTO>;
	       while( $lectura_linia !~ m/ALIGNMENTS/ ){
	            $lectura_linia = <ALINEAMIENTO>;
	            if( $lectura_linia !~ m/ALIGNMENTS/ ){
	                @linia_camps = ( $lectura_linia =~ m/\w+/g );  #llegir la linia per camps
	                $alineament[$linia_camps[0]][0]=$linia_camps[1];
	            }
	       }
	    }


	    ## LLEGEIX LA SEQ�NCIA DE CADA POSICIO I L'ASSIGNA A LA SEVA PROTEINA
            
	    if ($lectura_linia =~ m/ALIGNMENTS/ ){
	        while ($lectura_linia =~ m/ALIGNMENTS/ ){
	            @interval_prots = ( $lectura_linia =~m/\d+/g );
	            #print "minimo: $interval_prots[0]\n ";
	            #print "maximo: $interval_prots[1]\n ";
                    undef @aa;undef @ss;undef @bb;undef @s;undef @b;
	            $lectura_linia = <ALINEAMIENTO>;
                    if ($lectura_linia =~ m/ALIGNMENTS/){undef @aa;undef @ss;undef @bb;undef @s;undef @b;}
	            while( ($lectura_linia !~ m/ALIGNMENT/) && ($lectura_linia !~ m/SEQUENCE\s+PROFILE/) ){

	                $lectura_linia = <ALINEAMIENTO>;
                        if ($lectura_linia =~ m/ALIGNMENTS/){undef @aa;undef @ss;undef @bb;undef @s;undef @b;}
 	                if( ($lectura_linia !~ m/ALIGNMENT/) && ($lectura_linia !~ m/SEQUENCE\s+PROFILE/) ){
                              $step++;
                               $aa=substr($lectura_linia,14,1);
                               $s=substr($lectura_linia,17,1); 
                               if ($s eq " " || $s eq "-"){
                                push @ss,"C";
                               }else{
                                push @ss,$s; 
                               }
                               $b=substr($lectura_linia,36,3);
                               if ($sol->{$aa}){
                                 $bb=int(10*$b/$sol->{$aa});
                                 if ($bb < 0){$bb=0;}
                                 if ($bb >= 10) {$bb= 9;}
                               }
                               push @bb,$bb;


                            
	                    #####!!!!!!!!!!!Aquestes dos seguents linies es poden posar en una sola!!!he de consultar el llibre
                            
	                    @linia_camps = ( $lectura_linia =~ m/\w+/g );  #llegir la linia per camps
	                    $posicio = $linia_camps[0];    #identifica la posicio que s'esta llegint en aquell moment
	                    @linia_camps = ( $lectura_linia =~ m/[A-Z|a-z\s+.]+(?=\n)/g );  #selecciona de la linia aquells fragments que es corresponen a sequencia


	                    if ( $linia_camps[0] ){        #nom� entra si en aquella linia hi ha sequencia
	                        @secuencia = split( //, $linia_camps[0]);   #separa la sequencia per caracters

	                         #assigna els caracters a la seva posicio adequada a la matriu final de l'alineament
	                         for($i=$interval_prots[0]; $i<=$interval_prots[1]; $i++){
	                                $alineament[$i][$posicio]=$secuencia[$i-$interval_prots[0]+2];
	                         }
	                    }
	                }
	            }
	        }

	    }
	}

	close ALINEAMIENTO;

    $alineament[0][0]="Proteines";

	##OMPLE ELS BUITS AMB PUNTS ALL�ON CALGUI

	for ($i=1; $i<scalar(@alineament); $i++){
	    for ($j=1; $j<=$posicio; $j++){
	        if( !$alineament[$i][$j] || $alineament[$i][$j]eq" " ){
	                $alineament[$i][$j]=".";
	            }
	    }
	}

	@alineament = separa_matrius_HSSP(\@alineament);

	##treiem els limits inutils dels alineaments

        undef @limit_inf;
        undef @limit_sup;

	for( $i=0; $i<scalar(@alineament); $i++ ){
                printf "Alignment I $i\n";
		@{$alineament[$i]} = treure_limits_HSSP(\@{$alineament[$i]});
	}

	return @alineament;





		##FUNCIO QUE SEPARA UN ALINEAMENT LLEGIT DIRECTAMENT EN UN FITXER HSSP EN TOTS ELS SUBALINEAMENTS QUE CONT�		##RETORNA UN VECTOR ON A CADA POSICI�CONT�UNA MATRIU AMB CADASCUN DELS ALINEAMENTS

		sub separa_matrius_HSSP{

		my @origen = @{$_[0]};

		my $i;
		my $j;
		my $k;

		my @temporal;
		my @blocs;
		my @alineaments;

		my $inicio=0;
		my $fin=0;

		my $trobat;

		for( $i=1; $i<scalar(@origen); $i++){
			for( $j=1, $inicio=0, $fin=0; $j<scalar(@{$origen[$i]}); $j++ ){
			if( !$inicio ){
				if( !($origen[$i][$j] eq ".") ){
				$inicio = $j;
				}
			}
			else{
				if( !($origen[$i][$j] eq ".") ){
				$fin = $j;
				}
			}
			}
			$temporal[$i][0]=$inicio;
			$temporal[$i][1]=$fin;
		}

		#IMPRIMEIX PER PANTALLA ELS LIMITS
		#for( $i=1; $i<scalar(@temporal); $i++){
		#	for( $j=0; $j<2; $j++){
		#	    print "$temporal[$i][$j]  ";
		#	}
		#	print "\n";
		#    }

		#DETECTA BLOCS DE LIMITS
		$blocs[0][0]=$temporal[1][0];
		$blocs[0][1]=$temporal[1][1];

		for( $i=1; $i<scalar(@temporal); $i++){
			for( $j=0; $j<2; $j++){
			for( $k=0, $trobat=0; $k<scalar(@blocs) && !$trobat ; $k++){

				if( $temporal[$i][0] >= $blocs[$k][0] && $temporal[$i][1] <= $blocs[$k][1] ){
				$trobat=1;
				#print "l'hem trobat tot a dins!\n";
				}
				elsif( $temporal[$i][0] < $blocs[$k][0] &&  $temporal[$i][1] >= $blocs[$k][0] ) {
				$blocs[$k][0] = $temporal[$i][0];
				$trobat=1;
				#print "trobat!\n";
				if( $temporal[$i][1] > $blocs[$k][1] ){
					$blocs[$k][1] = $temporal[$i][1];
				}
				}
				elsif( $temporal[$i][1] > $blocs[$k][1] && $temporal[$i][0] <= $blocs[$k][1] ) {
				$blocs[$k][1] = $temporal[$i][1];
				$trobat=1;
				#print "trobat en 2\n";
				if( $temporal[$i][0] < $blocs[$k][0] ){
					$blocs[$k][0] = $temporal[$i][0];
				}
				}

			}

			if( !$trobat ){
				$blocs[$k][0]=$temporal[$i][0];
				$blocs[$k][1]=$temporal[$i][1];
			}
			}
		}


		#IMPRIMEIX PER PANTALLA ELS BLOCS
		#    for( $i=0; $i<scalar(@blocs); $i++){
		#	for( $j=0; $j<2; $j++){
		#	    print "$blocs[$i][$j]  ";
		#	}
		#	print "\n";
		#    }


		my $p=0;

		for( $i=1; $i<scalar(@temporal); $i++){
			for( $k=0, $trobat=0; !$trobat ; $k++){
			if( $temporal[$i][0] >= $blocs[$k][0] && $temporal[$i][1] <= $blocs[$k][1] ){
				#print "l'hem trobat a dins!!!!\n";
				$trobat=1;
				if( $alineaments[$k] ){
				$p=scalar(@{$alineaments[$k]});
				}
				else{
				$p=0;
				}
				@{$alineaments[$k][$p]} = @{$origen[$i]};
				$p++;
			}
			}

		}

		return @alineaments;

		}



		##FUNCIO PER TREURE ELS FRAGMENTS LIMITS SENSE SEQUENCIA EN ELS ALINEAMENTS SEPARATS OBTINGUTS DELS HSSP
		sub treure_limits_HSSP {

			my @alineament = @{$_[0]};
			my @alignretorn;

			my $i; my $j; my $k;

			my $surt;

			my $limit_inf=1000000;
			my $limit_sup = 0;

			##detecta el limit inferior (SEMPRE SER�UNA UNITAT SUPERIOR A LA REAL!!!)
			for( $i=0; $i<scalar(@alineament); $i++ ){
				for( $j=1, $surt=0; $j<scalar(@{$alineament[$i]}) && !$surt; $j++ ){
					if( !($alineament[$i][$j] eq "." ) ){
						$surt=1;
					}
				}
				if( $j<$limit_inf ){
					$limit_inf = $j;
				}
			}

			##detecta el limit superior
			for( $i=0; $i<scalar(@alineament); $i++ ){
				for( $j=1; $j<scalar(@{$alineament[$i]}); $j++ ){
					if( !($alineament[$i][$j] eq ".") ){
						$surt = $j;
					}
				}
				if( $surt>$limit_sup ){
					$limit_sup = $surt;
				}
			}

			##reescriu l'alineament 

                        push @limit_inf,$limit_inf-1;
                        push @limit_sup,$limit_sup;
			for($i=0; $i<scalar(@alineament); $i++ ){
				$alignretorn[$i][0] = $alineament[$i][0];
				for( $j=$limit_inf-1, $k=1; $j<=$limit_sup; $j++, $k++ ){
					$alignretorn[$i][$k] = $alineament[$i][$j] ;
				}
			}
			return @alignretorn;
		}
}



##LLEGEIX MATRIU de l'alineament i l'escriu en un fitxer de sortida
## Par�etres d'entrada: 
##	0: la matriu de l'alineament; 
##	1: el nom del fitxer de sortida
sub write_align{

	my @alineament = @{$_[0]};
	my $i;
	my $j;
	my $fitxer_sortida = $_[1];

	open(NEWALIGNMENT, ">$fitxer_sortida" ) || die "no es pot obrir el fitxer: $!";
	for ($i=0; $i<scalar(@alineament); $i++){
	    print NEWALIGNMENT "$alineament[$i][0]\t";
	    if( length($alineament[$i][0])<8 ){ print NEWALIGNMENT "\t"; }  #arregla el problema dels tabuladors en noms de proteina petits
	        for ($j=1; $j<scalar(@{$alineament[$i]}); $j++){
	            if( $alineament[$i][$j] ){
	                print NEWALIGNMENT "$alineament[$i][$j]";
	            }
	        }
	    print NEWALIGNMENT "\n";
	}

	close NEWALIGNMENT;
}

##LLEGEIX MATRIU de l'alineament i l'escriu en un fitxer de sortida
## Par�etres d'entrada: 0: la matriu de l'alineament; 1: el nom del fitxer de sortida; 2: quin: com cada fila � un alineament, indiquem quin d'aquests volem imprimir
sub write_align_HSSP{

	my @alineament = @{$_[0]};
	my $i;
	my $j;
	my $fitxer_sortida = $_[1];
	my $quin = $_[2];

	open(NEWALIGNMENT, ">$fitxer_sortida" ) || die "no es pot obrir el fitxer: $!";
	for ($i=0; $i<scalar(@{$alineament[$quin]}); $i++){
	    print NEWALIGNMENT ">M1 $alineament[$quin][$i][0]\n";
            undef @aa;
            undef @sp;
            undef @bp;
            undef @sm;
            undef @bm;
	        for ($j=1; $j<scalar(@{$alineament[$quin][$i]}); $j++){
	            if( $alineament[$quin][$i][$j] ){
                        if ($alineament[$quin][$i][$j] eq "."){
                            push @aa,"-";push @sp,"-";push @bp,"-";
                        }else{

                            push @aa,uc($alineament[$quin][$i][$j]);
                            push @sp,$ss[$j+$limit_inf[$quin]-2];
                            push @bp,$bb[$j+$limit_inf[$quin]-2];
                        }
                        push @sm,$ss[$j+$limit_inf[$quin]-2];
                        push @bm,$bb[$j+$limit_inf[$quin]-2];
	            }
	        }
              $n=0;
              for $i (0..$#aa){
               print NEWALIGNMENT "$aa[$i]";
               $n++;
               if ($n >79 ){print NEWALIGNMENT "\n";$n=0;}
              }
	    print NEWALIGNMENT "\n";


 	    $data={};
  	    $data={   AA  => [@aa],
       	      	      SS  => [@sm],
        	      ACC => [@bm],};

             undef @mizu;
 	     (@mizu) = &Mizuguchi($data);

 	     if (defined $matrix_s){
  	     printf NEWALIGNMENT "#%8.5f  %s\n",$weight,$matrix_s;
   	    $n=0;
   	    for $i (0..$#sp) {
   	     print NEWALIGNMENT "$sp[$i]";
   	     $n++;
   	     if ($n >79 ){print NEWALIGNMENT "\n";$n=0;}
   	    }
  	     print NEWALIGNMENT "\n";
  	    }
  	    if (defined $matrix_b){
  	     printf NEWALIGNMENT "#%8.5f  %s\n",$weight,$matrix_b;
  	     $n=0;
  	     for $i (0..$#bp) {
   	     print NEWALIGNMENT "$bur->{$bp[$i]}";
   	     $n++;
    	    if ($n >79 ){print NEWALIGNMENT "\n";$n=0;}
   	    }
  	     print NEWALIGNMENT "\n";
  	    }
  	    if (defined $matrix_m){
  	     printf NEWALIGNMENT "#%8.5f  %s\n",$weight,$matrix_m;
   	    $n=0;
   	    for $i (0..$#mizu) {
   	     if ($sp[$i] eq "-") {print NEWALIGNMENT "-";}else{print NEWALIGNMENT "$mizu[$i]";}
    	    $n++;
   	     if ($n >79 ){print NEWALIGNMENT "\n";$n=0;}
   	    }
   	    print NEWALIGNMENT "\n";
  	    }

	}

	close NEWALIGNMENT;
}

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
    if ($i>0){$ii=$i-1;}else{$ii=$i;}
    if ($strand > 0 ) { $beta++; }
    if ($helix  > 0 ) { $alfa++; }
    if ($beta == 1) {$b0=$ii;$acc_beta=$bb[$ii-1];}
    if ($alfa == 1) {$a0=$ii;if ($bb[$ii-1]>5) {$acc_alfa=1;}else{$acc_alfa=0;} }
    if ($strand > 0 && $beta > 1){ $acc_beta+=$bb[$ii-1];   }    
    if ($helix > 0 && $alfa > 1){ if ($bb[$ii-1]>5) {$acc_alfa+=1;} }
    if ($strand == 0 && $beta > 0){
        $size=$ii-$b0;
        if ($size < 1) {
         $beta=0;
         push @mizu,"C";
         if ($ss[$ii] ne "G" && $ss[$ii] ne "H" && $ss[$ii] ne "E" && $ss[$ii] ne "") {push @mizu,"C";}
        }else{
         $acc_beta=$acc_beta/$size;
         $short=$medium=$large=0;
         if ($size <= 3){$short=1;}
         if ($size <= 8 && $size > 3){$medium=1;}
         if ($size > 8) {$large=1;}
         if ($short == 1){if ($acc_beta <=2){$type="D";}else{$type="S";}}
         if ($medium== 1){if ($acc_beta <=2){$type="B";}else{$type="E";}}
         if ($large == 1){if ($acc_beta <=2){$type="M";}else{$type="N";}}
         for $j ($b0..($ii-1)){push @mizu,$type; }
         $beta=0;
         if ($ss[$ii] ne "G" && $ss[$ii] ne "H" && $ss[$ii] ne "E" && $ss[$ii] ne "") {push @mizu,"C";}
        }
       }
    if ($helix == 0 && $alfa > 0){
        $size=$ii-$a0;
        if ($size < 1) {
         $alfa=0;
         push @mizu,"C";
         if ($ss[$ii] ne "G" && $ss[$ii] ne "H" && $ss[$ii] ne "E" && $ss[$ii] ne "") {push @mizu,"C";}
        }else{
         $short=$large=0;
         if ($size <= 7){$short=1;}
         if ($size > 7) {$large=1;}
         if ($short == 1){if ($acc_alfa < 3){$type="A";}else{$type="F";}}
         if ($large == 1){if ($acc_alfa/2 < 3){$type="H";}else{$type="I";}}
         for $j ($a0..($ii-1)){push @mizu,$type;}
         $alfa=0;
         if ($ss[$ii] ne "G" && $ss[$ii] ne "H" && $ss[$ii] ne "E" && $ss[$ii] ne "") {push @mizu,"C";}
        }
       }
    if ($h310 > 0) {push @mizu,"G";}
    if ($ss[$i] eq "G") {$h310=1;$strand=0;$helix=0;}
    if ($ss[$i] eq "E") {$h310=0;$strand=1;$helix=0;}
    if ($ss[$i] eq "H") {$h310=0;$strand=0;$helix=1;}
    if ($ss[$i] ne "G" && $ss[$i] ne "H" && $ss[$i] ne "E" && $ss[$i] ne "") {
        $h310=0;$strand=0;$helix=0;
        if ($alfa==0 && $beta == 0) {push @mizu,"C";if ($i eq $#ss){$end=1;}}
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

