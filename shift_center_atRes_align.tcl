proc allign {comp ref molidmov} {
 	set transformation_matrix [measure fit $comp $ref]
	set movsel [atomselect $molidmov all]

	return  $transformation_matrix
}


set molidtraj [molinfo top get id]

set num_steps [molinfo $molidtraj get numframes]

set sel_residue [atomselect $molidtraj "protein and resid $cent_res and sidechain"]

set all [atomselect $molidtraj "all"]



set protein_i [atomselect $molidtraj "protein"]


for {set frame 0} {$frame < $num_steps} {incr frame} {
#centering
		$sel_residue frame $frame
		$all frame $frame
		set shift_vec1 [measure center $sel_residue weight mass]

		$all moveby [vecinvert $shift_vec1]
		set shift_vec1 [measure center $sel_residue weight mass]
		$all moveby [vecinvert $shift_vec1]

#aligning
		#$all frame $frame
		#$sel_residue frame $frame
		#set aligner [atomselect $molidtraj "protein and resid $cent_res and sidechain" frame first]
		#set aligner [atomselect $molidtraj "protein" frame 0]
		#$protein_i frame $frame 
		#set M [measure fit $protein_i $aligner]	 
		#$all move $M	
	
}




