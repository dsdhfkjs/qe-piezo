#!/bin/bash
#SBATCH -p A
#Script for either creating starting directories of a material or visualizing the relaxed unit cell

material_arr=("BN" "AlN" "GaN" "InN")
piezoversion_arr=("5.4-smart" "6.0-smart" "6.0-smart-kgrid")
dosversion_arr=("5.4-smart-k812" "6.0-smart-k812" "6.0-smart-kgrid") #
for material in ${material_arr[@]}; do

#System settings and conversion factors
filename=$(echo $material | awk '{print tolower($1)}')
convpath="Piezo/6.0-smart/"
convfile="$filename.vcrelax.out"

if [[ ! -d $material || ! -d $material/Standard || ! -d $material/Bands || ( ! -z $dosversion_arr && ! -d $material/DOS ) ]]; then #Create directories

	if [[ ! -d $material ]]; then mkdir $material; fi

	if [[ ! -d $material/Standard ]]; then mkdir $material/Standard; fi

	if [[ ! -d $material/Piezo || ! -z $piezoversion_arr ]]; then
		if [[ ! -d $material/Piezo ]]; then mkdir $material/Piezo; fi
	fi

	if [[ ! -z $piezoversion_arr ]]; then
		for piezoversion in ${piezoversion_arr[@]}; do
			if [[ ! -d $material/Piezo/$piezoversion ]]; then mkdir $material/Piezo/$piezoversion; fi
		done
	fi

	if [[ ! -z $dosversion_arr && ! -d $material/DOS ]]; then
		mkdir $material/DOS
		for dosversion in ${dosversion_arr[@]}; do
			if [[ ! -d $material/DOS/$dosversion ]]; then mkdir $material/DOS/$dosversion; fi
		done
	fi

	echo "Created directories for $material."

elif [[ -a "$material/Standard/$convfile" || -a "$material/$convpath$convfile" ]]; then #Generate visualizations of the unit cell
		
	#Run settings
	module load gcc/4.9.2
	
	#Save directory
	cd "./$material/Standard"
		
	#Ready the reference file
	if [[ ! -a $convfile && -a "../$convpath$convfile" ]]; then cp -p "../$convpath$convfile" $convfile ; fi
	
	#Open the PWSCF output file in XCrySDen then save as $filename.xsf
	xcrysden --pwo $convfile
	
	#Open .xsf in VESTA and save as .vesta and take a screenshot by exporting the vector image
	VESTA $filename.xsf
	
	cd ../../

fi #for making directories or visualizations
done #per material
