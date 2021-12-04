proc writedat1 {x y filename} {
	set outfile [open "$filename.txt" w]
	puts $outfile "frame \t RMSD"
	for {set i  0} {$i < [llength $x]} {incr i} {
		
		puts $outfile "[lindex $x $i] \t [lindex $y $i]"
	}
        
	close $outfile
}


proc writedat2 {x y filename} {
	set outfile [open "$filename.dat" w]
	
	puts $outfile "frame \t Y1 \t Y2 \t Y3"
		for {set i  0} {$i < [llength $x]} {incr i} {
			set yi [lindex $y $i]
			puts $outfile "[lindex $x $i] \t [lindex $yi 0] \t [lindex $yi 1] \t [lindex $yi 2]"
		}
	
			
        close $outfile
}

proc frame_to_time {frame md_step sampling_rate} { 
#md_step in fs 
	set t [expr [expr $md_step*$frame*$sampling_rate]/1000 ]
	return $t
 
#output in ns

}

#set path /home/sk87/Programs/Expanse/Expanse_output/ctpr8_free_post500_res35
#set filename ctpr8_wb_eq_35_exc





proc COMvst {state DCD PSF Ns Nf resid centres}	{
	set filename ctpr8_${state}_${resid}
	set mol [mol new $PSF type psf waitfor all]	
	mol addfile $DCD type dcd first $Ns last $Nf step 1 waitfor all
		
		
 

	set num_frames [molinfo top get numframes]

	set cent_res $centres
	


	source /home/sk87/Programs/DNAPolymerase2021/Phi29_CTPR8_PDBs/analysisscritps/movie_making/shift_center_atRes.tcl

	set num_frames [molinfo top get numframes]
	
	set first_frame $Ns
	set last_frame  $Nf
	set stride       2500

	 
	set sel_resID [atomselect top "protein and resid $resid and sidechain"]
	set atom_names [$sel_resID get name]
	
	#sampling rate( in ns) 
	set steps_per_frame 2500  
	#stpesize in femtoseconds
	set time_per_step   4.0    

	set sampling_time   [expr [expr $steps_per_frame*$time_per_step]/1000] 
	#in ps 
	set avg_time        $sampling_time  
	
	set steps_for_avg   [expr int([expr [expr 1.0/$steps_per_frame] *[expr [expr $avg_time *1000]/$time_per_step]])]
		
			#calculating start and end time
	set start_time [expr int([frame_to_time $first_frame $time_per_step $steps_per_frame])]
	set end_time   [expr int([frame_to_time $last_frame $time_per_step $steps_per_frame])]
	
	set outdir  "COM_Q_${state}_rid${resid}_${avg_time}ns_${start_time}_${end_time}"
		
	file mkdir $outdir
	cd $outdir
	#foreach $outdir [glob *] {
    	#	file delete -force -- $path
	#}
	
	set Q {}
	
	foreach name $atom_names {
			set POS {}
			set t {}

			
			
			for {set i 0} {$i < $num_frames} {incr i $steps_for_avg} {
				set seli [atomselect top "protein and resid $resid and sidechain and name $name" frame $i ]
				set posi { 0 0 0 }
				for {set j 0} { $j < $steps_for_avg} {incr j} {
					set posj [$seli get {x y z}]
					
					set posj [lindex $posj 0]
					
					set posi [vecadd $posi [vecscale $posj [expr 1.0/$steps_for_avg]]]	
			
				}
				
				set $POS  [lappend POS   $posi ]
				set $t    [lappend t [expr $i+$Ns ]]
				
				
			}
			
			set q [ [atomselect top "protein and resid $resid and sidechain and name $name"] get charge ] 
			
			set $Q [lappend Q $q]


			
			set POS_file      Pos_rid${resid}_${name}


			writedat2 $t $POS ${POS_file} 
			
			
	}
	writedat1 $atom_names $Q Charge
	set GEOM_out [pwd]
	
	
}

