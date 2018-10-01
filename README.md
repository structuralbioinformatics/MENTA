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


