#!/bin/bash
CONFIG=${1? Error: no configuration file entered}
echo $CONFIG
chmod u+x $CONFIG
source $CONFIG

CD=$(echo $PWD)
mkdir $OUTFILE
rm -r "$OUTFILE"/*
cd $OUTFILE
>error.log
DCD_LEN=$(catdcd -num $DCDFILE | grep -i "Total Frames" | tr -d -c 0-9)

CHUNK=5000

##$((DCD_LEN))
for (( i=0; i<$((DCD_LEN)); i+=5000 ))
do
	Ns=$i
	if [ $((i+CHUNK)) -gt $((DCD_LEN)) ] 
	then
		Nf=-1
	else
		Nf=$((i+CHUNK))
	fi
	
	find $V_DIR -type d -name "*_$Ns_$Nf" | head -n 1 
	V_OUT=$(find $V_DIR -type d -name "*_$Ns_$Nf" | head -n 1)
	echo $V_OUT
	if [ $STATE == "Charged" ]
	then
	
		printf "$STATE\n$Ns\n$Nf\n$DCDFILE\n$PSFFILE\n$CHARGEDRES\n$CHARGEDRES"| vmd -dispdev text -e run_script.tcl | tee temp.out
	else
		printf "$STATE\n$Ns\n$Nf\n$DCDFILE\n$PSFFILE\n$NEUTRES\n$CHARGEDRES"| vmd -dispdev text -e /home/sk87/Programs/DNAPolymerase2021/Phi29_CTPR8_PDBs/analysisscritps/Residue_geom/run_script.tcl | tee temp.out
	fi
	
        GEOM_OUT=$(grep "COM_Q" temp.out)
	echo $GEOM_OUT
	rm temp.out
	
	function useminiconda() { PATH="/home/cmaffeo2/miniconda3/bin:$PATH";};
	useminiconda
	export PATH
	python3 /home/sk87/Programs/DNAPolymerase2021/Phi29_CTPR8_PDBs/analysisscritps/PMEpot/energy_gap_parallel.py $STATE $V_OUT $GEOM_OUT $RES_TYPE $CHARGEDRES $NEUTRES $FILENAME 2>> error.log
done
mkdir X_"$OUTFILE"
mv *.dat X_"$OUTFILE"
cd $CD
