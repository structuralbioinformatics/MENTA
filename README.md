# MENTA

Program for the alignment of profiles, creation of profiles and multiple sequence alginments using sequences
Folder perl contains scripts to help the creation and use of inputs using secondary structure from DSSP or predicted by SABLE
or profiles from HSSP
Format of the input is FastA like. If the title starts by >M#, with # integer, the sequences are grouped into a profile defined as M#

PARAMETERS
    
	 -i 	 File with Sequences      
	 -o 	 Output File     
	 -op	 Minimum number of sequences on profiles to write in the output     
	 -a 	 Elements (the first two for gap/extension)<default:--ABCDEFGHIKLMNPQRSTVWXYZ>      
	 -w 	 File with Weights and Substitution-Matrix files      
	 -wj	 File with Weights for Scoring Method (1..13 weights)     
	 -n 	 Minimum Alignment <default: 2>     
	 -r 	 Tolerance Error <default: 1.0e-6>     
	 -e 	 Tolerance Extension  <default: 10 positions>     
	 -ej	 Tolerance Extension for local alignment     
	    	 only if Scoring-Method "0" is applied  <default: NULL => 500 positions >     
	 -ep	 Tolerance Extension in score penalty  <default: 100.0 >     
	 -ejp	 Tolerance Extension penalty for local alignment     
	    	 only if Scoring-Method "0" is applied  <default: NULL => 500*median-score>     
	 -f 	 Tolerance Factor <default: 0.0>     
	 -s 	 Tolerance Score <default: -100.0>     
	 -id	 Minimum Ratio of sequence identity to create a new profile     
	 -homo	 Minimum Ratio of sequence similarity to create a new profile     
	 -cluster	 Maximum number of profiles with common sequences     
	 -gid	 Gradient factor to apply on restrictions of ID ratio in each iteration    
	 -ghom	 Gradient factor to apply on restrictions of HOMO ratio in each iteration    
	 -psic	 Pseudo Counting flag: ON (simplified approach) / OFF (full pseudo-counting) <default OFF>     
	 -sub	 Level of Global subalignments <default: 3>     
	    	   0. Get all subalignments      
	    	   1. Switch off additional pairs of gaps in the same position     
	    	   2. Switch off additional gaps in one direction     
	    	   3. Switch off all additional gaps (discard subalignments)     
	 -j 	 Scoring Method <default: 0>     
	    	   0. MEntA     
	    	   1. Sum of pairs with Raw Frequencies     
	    	   2. Sum of pairs with Efficient Frequencies (F)      
	    	   3. Sum of pairs with Target Frequencies (Q)      
	    	   4. Maximum pair      
	    	   5. Pearson Correlation of Efficient Frequencies (F)     
	    	   6. Pearson Correlation of Target Frequencies (Q)     
	    	   7. Pearson Correlation of Natural logarithm of F/q     
	    	   8. Pearson Correlation of Natural logarithm of Q/q     
	    	   9. PICASSO     
	    	  10. PICASSO-Q     
	    	  11. COMPASS (unweighted)     
	    	  12. COMPASS (weighted)     
	    	  13. PROFSIM      
	 -fmt 	 Output Format <default: 0>     
	    	   0. MEntA format     
	    	   1. PIR format       
	    	   2. CLUSTAL format      
	    	  10. MEntA format (includes the alignment of properties)     
	    	  11. PIR format (includes the alignment of properties)     
	    	  12. CLUSTAL format (includes the alignment of properties)     
	 -evd 	 Method of E-value calculation <default: 0>     
	    	   0. MEntA      
	    	   1. Island     
	    	   2. Murzin P-value     
	    	   3. Dummy     
	 -g 	 Global Mode      
	 -l 	 Local Mode      
	 -h 	 Print help 

============================================================

EXAMPLE:

GO TO "example" folder, where it contains the sequence and profile files required to run the example

1) Align and create local motif profiles using single sequences

We assume you have installed SABLE and the program runs as run.sable.
The output of SABLE is "input_fasa"+SABLE.out (i.e run.sable src8.fa produces src8.faSABLE.out)

SABLE/run.sable src8.fa  
perl ../perl/SABLE2Menta.pl -i src8.faSABLE.out -f src8 -m ../matrix/Mizuguchi.dat  -o src8
SABLE/run.sable ubx8.fa
perl ../perl/SABLE2Menta.pl -i ubx8.faSABLE.out -f ubx8 -m ../matrix/Mizuguchi.dat  -o ubx8
SABLE/run.sable faf1.fa
perl ../perl/SABLE2Menta.pl -i faf1.faSABLE.out -f faf1 -m ../matrix/Mizuguchi.dat  -o faf1

cat src8.menta > iMotif.menta
cat ubx8.menta >> iMotif.menta
cat faf1.menta >> iMotif.menta

echo " 0.80          ../matrix/blosumG.dat" > matrix.data
echo " 0.10          ../matrix/blosum.dat" >> matrix.data
echo " 0.10          ../matrix/pam.dat"    >> matrix.data

echo " 14.0" > input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat
echo "  0.0" >> input_method_WJ.dat
echo " 14.0" >> input_method_WJ.dat


../bin/menta -i iMotif.menta -w matrix.data -wj input_method_WJ.dat -n 10  -fmt 10 -j 0 -s -100.0 -id 0.3 -homo 0.9 -gid 0.5 -ghom 0.9 -cluster 20 -l -evd 1 -op 4  -o test_iMotif.out > test_iMotif.log

sed -e "s/0.00000/0.25000/g" iMotif.menta > iMotif_ss.menta

../bin/menta -i iMotif_ss.menta -w matrix.data -wj input_method_WJ.dat -n 10  -fmt 10 -j 0 -s -100.0 -id 0.3 -homo 0.9 -gid 0.5 -ghom 0.9 -cluster 20 -l -evd 1 -op 4  -o test_iMotif_ss.out > test_iMotif_ss.log



2) Use HSSP files to create two profiles and compare them

perl ../perl/HSSP2MEntA.pl -i 1aw0.hssp -m ../matrix/Mizuguchi.dat -w 0.5 -o 1aw0
perl ../perl/HSSP2MEntA.pl -i 1kdc.hssp -m ../matrix/Mizuguchi.dat -w 0.5 -o 1kdc

cat 1aw0_0.menta > test_1aw0_1kdc.menta
sed -e "s/>M1/>M2/g" 1kdc_0.menta >> test_1aw0_1kdc.menta

../bin/menta -i test_1aw0_1kdc.menta  -w matrix.data -wj input_method_WJ.dat -n 10  -fmt 10 -j 0 -s -100.0 -id 0.3 -homo 0.9 -gid 0.5 -ghom 0.9 -cluster 20 -g -evd 1 -op 4  -o test_HSSP.out >  test_HSSP.log


