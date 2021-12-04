The code requires python and vmd on the system which is running. 

Specify the path to the psf, dcd and folder iwth .dx file for the potential energy landscape for every frame in the energy_gap.config file. 

Run the script using bash commond ./energygap_workflow.sh energy_gap.config. This integrates the vmd tcl scripts which generate the postiion data for the residue and the python script which calculates the energy gap as a function of time. 
The results are stored in the .dat file inside the output file specified in the configuration file 
