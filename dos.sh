#!/bin/bash
#SBATCH -p A

#Initial directory
cd BaNiO3/DOS/6.0-smart

#Run settings
path=$PATH
export OMP_NUM_THREADS=1
runprefix="mpirun -np 4 pw.x" #salloc -p A -n 4

#System settings and conversion factors
pseudodir="/home/sand/pang/Research/Pseudo/"
forc=9.72344907103d-5
ewfc=73.49861765 #Ry
kptden=2000 #k-point density per reciprocal atom
kptden_dos=7000 #k-point density per reciprocal atom
Avogadro=6.022140857e23 #atoms/mol
BohrtoAng=0.529177249 #Angstrom/bohr
Ang3tocm3=1.0e-24 #cm^3/Angstrom^3
Rybohr3toGPa=14710.5164168 #GPa/(Ry/bohr^3)
RytoeV=13.605698066 #eV/Ry
eVtoJ=1.60217733e-19 #J/eV
boltz=1.38064852e-23 #J/K

#Select the materials to study
material_arr=("PbTiO3_c") # "AlN" "LaOF" "VFeSb" "YOF_c" "NiAsS" "BaNiO3" "BaTiO3_t"
for material in ${material_arr[@]}; do

if [[ $material = 'BaNiO3' ]]; then
loop_arr=('6.0-kgrid' '6.0-kgrid-p') #'6.0-k812-p'
plotrange="-40to10"
overflow=0.75
crys_sys="hex"
point_grp="6mm"
convpath="../../Piezo/6.0/"
kptpathden=5
kptwt=10
kpthighsym="
0.000000 0.000000 0.000000 10
0.500000 0.000000 0.000000 10
0.333333 0.333333 0.000000 10
0.000000 0.000000 0.000000 10
0.000000 0.000000 0.500000 10
0.500000 0.000000 0.500000 10
0.333333 0.333333 0.500000 10
0.000000 0.000000 0.500000 10
0.500000 0.000000 0.500000 start
0.500000 0.000000 0.000000 end
0.333333 0.333333 0.000000 start
0.333333 0.333333 0.500000 end"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ni   58.6934	Ni.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'VFeSb' ]]; then
loop_arr=('6.0-smart-k812a' '6.0-smart-k812a-p')
plotrange="-40to10"
overflow=0.75
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="-43m"
kptpathden=5
kptwt=10 #GM - X - W - K - GM - L - U - W - L - K|U - X
kpthighsym="
0.000000 0.000000 0.000000 10
0.000000 0.500000 0.000000 10
0.250000 0.500000 0.000000 10
0.375000 0.375000 0.000000 10
0.000000 0.000000 0.000000 10
0.250000 0.250000 0.250000 10

0.?00000 0.?00000 0.?00000 start
0.?00000 0.?00000 0.?00000 end"
atomic_species="
V    50.9415	V.GGA-PBE-paw.UPF  
Fe   55.845	Fe.GGA-PBE-paw.UPF 
Sb   121.76	Sb.GGA-PBE-paw.UPF "

elif [[ $material = 'YOF_c' ]]; then
loop_arr=('6.0-smart-k812a' '6.0-smart-k812a-p')
plotrange="-40to10"
overflow=0.75
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="-43m"
atomic_species="
Y    88.90585	Y.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF 
F    18.998403	F.GGA-PBE-paw.UPF "

elif [[ $material = 'LaOF' ]]; then
loop_arr=('6.0-smart-k812' '6.0-smart-k812-p') #'5.4-smart' '6.0-smart-p'
plotrange="-36to10"
overflow=0.5
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="-43m"
atomic_species="
La   138.90547	La.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF 
F    18.998403	F.GGA-PBE-paw.UPF "

elif [[ $material = 'AlN' ]]; then
loop_arr=('6.0-smart-k812-p')
plotrange="-18to6"
overflow=0.5
convpath="../../Piezo/6.0-smart/"
crys_sys="hex"
point_grp="6mm"
atomic_species="
   Al	26.981539 	Al.GGA-PBE-paw.UPF  
   N	14.0067		N.GGA-PBE-paw.UPF "

elif [[ $material = 'BaTiO3_t' ]]; then
loop_arr=('6.0-smart-kgrid-p')
plotrange="-5to10,-40to10"
overflow=0.75
convpath="../../Piezo/6.0-smart/"
crys_sys="tet_high"
point_grp="4mm"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'KNbO3_c' ]]; then
loop_arr=('6.0-smart-k812' '6.0-smart-k812-p')
plotrange="-40to10"
point_grp="m-3m"
atomic_species="
K    39.0983	K.GGA-PBE-paw.UPF 
Nb   92.90638	Nb.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'PbTiO3_c' ]]; then
loop_arr=('6.0-smart-kgrid' '6.0-smart-kgrid-p') #6.0-smart-k812a
plotrange="-5to10,-40to10" #-10to10,
overflow=0.75
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="m-3m"
kptpathden=5
kptwt=10
kpthighsym="
0.000000 0.000000 0.000000 10
0.000000 0.500000 0.000000 10
0.500000 0.500000 0.000000 10
0.000000 0.000000 0.000000 10
0.500000 0.500000 0.500000 10
0.000000 0.500000 0.000000 10
0.500000 0.500000 0.000000 start
0.500000 0.500000 0.500000 end"
atomic_species="
Pb   207.2	Pb.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'GaPt' ]]; then
loop_arr=('6.0-smart-k812' '6.0-smart-k812-p')
plotrange="-40to10"
overflow=0.5
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="23"
atomic_species="
Ga   69.723	Ga.GGA-PBE-paw.UPF  
Pt   195.084	Pt.GGA-PBE-paw.UPF "

elif [[ $material = 'NaBrO3' ]]; then
loop_arr=('6.0-smart-k812' '6.0-smart-k812-p')
plotrange="-40to10"
overflow=0.5
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="23"
atomic_species="
Na   22.989769	Na.GGA-PBE-paw.UPF  
Br   79.904	Br.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'NiAsS' ]]; then
loop_arr=('6.0-smart-k812' '6.0-smart-k812-p')
plotrange="-40to10"
overflow=0.75
convpath="../../Piezo/6.0-smart/"
crys_sys="cubic"
point_grp="23"
atomic_species="
Ni   58.6934	Ni.GGA-PBE-paw.UPF  
As   74.9216	As.GGA-PBE-paw.UPF  
S    32.065	S.GGA-PBE-paw.UPF "
fi

#Finished runs

#VFeSb - Error in routine tweights (1): bad Fermi energy - how to solve
#YOF_c	- 
#NiAsS - Error in routine tweights (1): bad Fermi energy - how to solve

#Save directory and material name
cd ../../../$material
shdir=$(pwd); cd DOS/6.0-smart-k812/ 
filename=$(echo $material | awk '{print tolower($1)}')

for loop in ${loop_arr[@]}; do

#Load modules for the specific version
if [[ $loop = *"5.4-smart-k812"* ]]; then
	cd ../5.4-smart-k812/
	export PATH=/home/sand/pang/espresso-5.4.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/1.10.0
	kpt=8
	kptbig=12
elif [[ $loop = *"5.4-smart"* ]]; then
	cd ../5.4-smart/
	export PATH=/home/sand/pang/espresso-5.4.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/1.10.0
	kpt=6
	kptbig=10
elif [[ $loop = *"6.0-smart-kgrid"* ]]; then
	cd ../6.0-smart-kgrid/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=4
elif [[ $loop = *"6.0-smart-k812b"* ]]; then
	cd ../6.0-smart-k812b/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=12
elif [[ $loop = *"6.0-smart-k812a"* ]]; then
	cd ../6.0-smart-k812a/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=12
elif [[ $loop = *"6.0-smart-k812"* ]]; then
	cd ../6.0-smart-k812/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=12
elif [[ $loop = *"6.0-smart"* ]]; then
	cd ../6.0-smart/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=6
	kptbig=10
elif [[ $loop = *"6.0-kgrid"* ]]; then
	cd ../6.0-kgrid/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=4
elif [[ $loop = *"6.0-k812"* ]]; then
	cd ../6.0-k812/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptbig=12
elif [[ $loop = *"6.0"* ]]; then
	cd ../6.0/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=6
	kptbig=10
fi

#Check if we will plot
if [[ $loop = *"p"* ]]; then
	dos_plot="" #[all]
	pdos_plot="" #[all]
	bands_plot="[all]"
else
	dos_plot=""
	pdos_plot=""
	bands_plot=""
fi

#Set output file
if [[ $loop = *"p"* ]]; then output="$filename.$loop.$plotrange.txt"; else output="$filename.$loop.txt"; fi
echo -n "" > $output

#Material parameters
rundir=$(pwd); cd $shdir
cif="Standard/*conventional_standard.cif"
a_ref=$(grep 'cell_length_a' $cif | awk '{print $2}') #Angstrom
b_ref=$(grep 'cell_length_b' $cif | awk '{print $2}') #Angstrom
c_ref=$(grep 'cell_length_c' $cif | awk '{print $2}') #Angstrom
cosab_ref=$(python -c "from math import cos, radians; print cos(radians($(grep 'cell_angle_gamma' $cif | awk '{print $2}')))")
cosac_ref=$(python -c "from math import cos, radians; print cos(radians($(grep 'cell_angle_beta' $cif | awk '{print $2}')))")
cosbc_ref=$(python -c "from math import cos, radians; print cos(radians($(grep 'cell_angle_alpha' $cif | awk '{print $2}')))")
nat=$(sed -n '/_atom_site_occupancy/,$p' $cif | awk 'END {print NR-1}')
ntyp=$(grep 'chemical_formula_sum' $cif | awk '{print NF-1}')
pos_ref=$(tail -$nat $cif | awk '{print $1, $4, $5, $6}')
if [[ $loop = *"kgrid"* ]]; then
	read kpta_ref kptb_ref kptc_ref <<<$(python -c "import numpy as np; vec=np.array([$a_ref, $b_ref, $c_ref], float); print np.ceil(np.cbrt($kptden/($nat*((1/vec)[0])*((1/vec)[1])*((1/vec)[2])))*(1/vec))" | sed 's/.//;s/.$//;s/\.//g')
else
	kpta_ref=$kpt
	kptb_ref=$kpt
	kptc_ref=$kpt
fi
cd $rundir
echo "Using a k-point grid of $kpta_ref x $kptb_ref x $kptc_ref for vc-relax" >> $output

#Obtain the reference unit cell with our settings
echo -e "\nCalculating the relaxed parameters for the unstrained unit cell of $material."
cat > ${filename}.vcrelax.in << EOF
&CONTROL
		 calculation = 'vc-relax',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'from_scratch',
	       forc_conv_thr = $forc,
 /
 &SYSTEM
                       ibrav = 14,
                   	   A = $a_ref,
                   	   B = $b_ref,
                   	   C = $c_ref,
                       cosAB = $cosab_ref,
                       cosAC = $cosac_ref,
                       cosBC = $cosbc_ref,
                         nat = $nat,
                        ntyp = $ntyp,
                     ecutwfc = $ewfc,
                 occupations = 'fixed' ,
                     degauss = 0.00
 /
 &ELECTRONS
                    conv_thr = 1.0e-6,
                 mixing_beta = 0.3 ,
 /
 &IONS
		ion_dynamics = 'bfgs'
 /
 &CELL
	       cell_dynamics = 'bfgs'
 /
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_ref

K_POINTS automatic 
  $kpta_ref $kptb_ref $kptc_ref   0 0 0 
EOF
convfile="$filename.vcrelax.out"
if [[ ( ! -a $convfile || -z $(grep 'PWSCF  ' $convfile) ) && -a $convpath$convfile ]]; then cp -p $convpath$convfile $convfile ; elif [[ ! -a $convpath$convfile ]]; then $runprefix pw.x < $filename.vcrelax.in > $convfile; fi
#if [[ $loop != *"p"* ]]; then $runprefix pw.x < $filename.vcrelax.in > $convfile; fi

alat=$(grep -n "CELL_PARAMETERS" $convfile | tail -1 | cut -f1 -d:)
read alat a1_i a2_i a3_i b1_i b2_i b3_i c1_i c2_i c3_i <<<$(echo $(sed -n "$((alat)),$((alat+3))p" $convfile) | awk '{gsub(")",""); print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}')
a_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($a1_i**2) + abs($a2_i**2) + abs($a3_i**2))")
b_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($b1_i**2) + abs($b2_i**2) + abs($b3_i**2))")
c_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($c1_i**2) + abs($c2_i**2) + abs($c3_i**2))")
alpha_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($b1_i*$c1_i + $b2_i*$c2_i + $b3_i*$c3_i)/($b_i*$c_i)))")
beta_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($c1_i*$a1_i + $c2_i*$a2_i + $c3_i*$a3_i)/($c_i*$a_i)))")
gamma_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($a1_i*$b1_i + $a2_i*$b2_i + $a3_i*$b3_i)/($a_i*$b_i)))")
pos_vc=$(grep -n "ATOMIC_POSITIONS" $convfile | tail -1 | cut -f1 -d:)
pos_vc=$(sed -n "$((pos_vc+1)),$((pos_vc+nat))p" $convfile)
mass_i=0.0
for row in $(seq 0 $(echo $atomic_species | awk '{print (NF/3 -1)}')); do
	mass_i=$(python -c "print $mass_i + $(echo $atomic_species | awk '{print $(row*3 + 2)}' row=$row)*$(grep -o "$(echo $atomic_species | awk '{print $(row*3 + 1)}' row=$row)" <<< "$pos_vc" | wc -l)/$Avogadro")
done
volume_i=$(python -c "print ($BohrtoAng*$alat)**3*($a1_i*($b2_i*$c3_i - $b3_i*$c2_i) + $a2_i*($b3_i*$c1_i - $b1_i*$c3_i) + $a3_i*($b1_i*$c2_i - $b2_i*$c1_i))")
density_i=$(python -c "print $mass_i/($volume_i*$Ang3tocm3)")
printf "a_i: %0.4f Ang\t\tb_i: %0.4f Ang\t\tc_i: %0.4f Ang\n" $a_i $b_i $c_i >> $output
echo -e "alpha_i: $alpha_i degrees\t\tbeta_i: $beta_i degrees\t\tgamma_i: $gamma_i degrees" >> $output
printf "volume_i: %0.4f Ang^3\n" $volume_i >> $output
printf "density_i: %0.3f g/cm^3\n\n\n" $density_i >> $output
if [[ $loop = *"kgrid"* ]]; then
	read kpta kptb kptc <<<$(python -c "import numpy as np; vec=np.array([$a_i, $b_i, $c_i], float); print np.ceil(np.cbrt($kptden/($nat*((1/vec)[0])*((1/vec)[1])*((1/vec)[2])))*(1/vec))" | sed 's/.//;s/.$//;s/\.//g')
else
	kpta=$kpt
	kptb=$kpt
	kptc=$kpt
fi
echo "Using a k-point grid of $kpta x $kptb x $kptc for scf" >> $output

nelec=$(grep 'number of electrons' $convfile | awk '{print $NF}' | tail -1)
nbnd=$(printf "%.0f\n" $(python -c "print $overflow*$nelec"))

#Delete previous $filename.save and save_temp
if [[ $loop != *"p"* && -d "$filename.save" && ! -z "$(ls -A $filename.save)" ]]; then chmod -R a+w ./$filename.save; rm -r ./$filename.save/*; fi
if [[ $loop != *"p"* && -d "save_temp" && ! -z "$(ls -A save_temp)" ]]; then chmod -R a+w ./save_temp; rm -r ./save_temp/*; fi


echo "Running the self-consistent field calculation for the optimized cell."
cat > $filename.scf.in << EOF
&CONTROL
		 calculation = 'scf',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'from_scratch',
	       forc_conv_thr = $forc,
		   verbosity = 'high',
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
                         nat = $nat,
                        ntyp = $ntyp,
			nbnd = $nbnd, 
                     ecutwfc = $ewfc,
                 occupations = 'fixed' ,
                     degauss = 0.00 ,
 /
 &ELECTRONS
                    conv_thr = 1.0e-6,
                 mixing_beta = 0.3 ,
 /
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_vc

K_POINTS automatic 
  $kpta $kptb $kptc   0 0 0 

CELL_PARAMETERS alat
$a1_i $a2_i $a3_i
$b1_i $b2_i $b3_i
$c1_i $c2_i $c3_i
EOF
scffile="$filename.scf.out"
nscffile="$filename.nscf.out"
pwbandsfile="$filename.pwbands.out"
if [[ $loop != *"p"* && ( ! -a $scffile || -z $(grep 'PWSCF  ' $scffile) || ! -a $nscffile || -z $(grep 'PWSCF  ' $nscffile) || ! -a $pwbandsfile || -z $(grep 'PWSCF  ' $pwbandsfile) ) ]]; then $runprefix pw.x < $filename.scf.in > $scffile; fi

#Get the runtime, and total energy
scfetot=$(grep '!' $scffile | tail -1 | awk '{print $(NF-1)}')
scffermi=$(grep "Fermi energy\|highest occupied" $scffile | awk '{print $(NF-1)}')
scftime=$(grep "PWSCF" $scffile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')
echo -e "Fermi/highest occupied level: $scffermi eV" >> $output
echo -e "Total energy after scf: $scfetot Ry\n" >> $output

#Get k-point grid for nscf
if [[ $loop = *"kgrid"* ]]; then
	read kpta_big kptb_big kptc_big <<<$(python -c "import numpy as np; vec=np.array([$a_i, $b_i, $c_i], float); print np.ceil(np.cbrt($kptden_dos/($nat*((1/vec)[0])*((1/vec)[1])*((1/vec)[2])))*(1/vec))" | sed 's/.//;s/.$//;s/\.//g')
else
	kpta_big=$kptbig
	kptb_big=$kptbig
	kptc_big=$kptbig
fi
echo "Using a k-point grid of $kpta_big x $kptb_big x $kptc_big for nscf" >> $output

#Save $filename.save from SCF
if [[ $loop != *"p"* && ! -d "save_temp" ]]; then mkdir save_temp; fi
if [[ $loop != *"p"* ]]; then cp -r $filename.save ./save_temp; fi


echo "Running the non-self-consistent field calculation for the optimized cell."
cat > $filename.nscf.in << EOF
&CONTROL
		 calculation = 'nscf',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'from_scratch',
	       forc_conv_thr = $forc,
		   verbosity = 'high',
		  wf_collect = .true.
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
                         nat = $nat,
                        ntyp = $ntyp,
			nbnd = $nbnd, 
                     ecutwfc = $ewfc,
                 occupations = 'tetrahedra' ,
                     degauss = 0.00 ,
 /
 &ELECTRONS
                    conv_thr = 1.0e-6,
                 mixing_beta = 0.3 ,
 /
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_vc

K_POINTS automatic 
  $kpta_big $kptb_big $kptc_big   0 0 0 

CELL_PARAMETERS alat
$a1_i $a2_i $a3_i
$b1_i $b2_i $b3_i
$c1_i $c2_i $c3_i
EOF
#if [[ $loop != *"p"* && ( ! -a $nscffile || -z $(grep 'PWSCF  ' $nscffile) || ! -a $dosfile || -z $(grep 'DOS  ' $dosfile) ) ]]; then $runprefix pw.x < $filename.nscf.in > $nscffile; fi
if [[ $loop != *"p"* && ( ! -a $nscffile || -z $(grep 'PWSCF  ' $nscffile) || ! -a $pwbandsfile || -z $(grep 'PWSCF  ' $pwbandsfile) ) ]]; then $runprefix pw.x < $filename.nscf.in > $nscffile; fi

#Get the Fermi energy in eV and runtime
fermi=$(grep "Fermi energy\|highest occupied" $nscffile | awk '{print $(NF-1)}')
nscftime=$(grep "PWSCF" $nscffile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')
echo -e "Fermi/highest occupied level: $fermi eV\n" >> $output


echo "Calculating the total density of states after nscf."
cat > $filename.dos.in << EOF
&DOS
		prefix = '$filename',
	       fildos = '$filename.dos',
	        DeltaE = 0.1,
		ngauss = 0,
	       degauss = 0.0,
/
EOF
dosfile="$filename.dos.out"
#if [[ $loop != *"p"* && ( ! -a $dosfile || -z $(grep 'DOS  ' $dosfile) ) ]]; then $runprefix dos.x < $filename.dos.in > $dosfile; fi
#if [[ $loop != *"p"* && ( ! -a $dosfile || -z $(grep 'DOS  ' $dosfile) ) ]]; then dos.x < $filename.dos.in > $dosfile; fi
if [[ $loop != *"p"* ]]; then dos.x < $filename.dos.in > $dosfile; fi

#Get the runtime
dostime=$(grep "DOS" $dosfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')


echo "Calculating the projected density of states after nscf."
cat > $filename.pdos.in << EOF
&PROJWFC
		prefix = '$filename',
	       filpdos = '$filename.pdos',
	        DeltaE = 0.1,
		ngauss = 0,
	       degauss = 0.0,
/
EOF
pdosfile="$filename.pdos.out"
if [[ $loop != *"p"* ]]; then projwfc.x < $filename.pdos.in > $pdosfile; fi

#Get the runtime
pdostime=$(grep "PROJWFC" $pdosfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')

#Initialize the k-point path for the band structure
banddiv=$(($(grep -o "start" <<< $kpthighsym | wc -l) + 1))
kptlines=$(echo $kpthighsym | awk '{print (NF/4 -1)}')
kpttot=$(( kptlines*kptpathden + 1 - (banddiv-1)*(kptpathden-1) ))
#echo -e "kpttot: $kpttot\tkptlines: $kptlines\tkptpathden: $kptpathden\tbanddiv: $banddiv"
kptpath="
$kpttot
$(echo $kpthighsym | awk '{print $1, $2, $3, $4}')"

#Generate the k-point path for the band structure
for row in $(seq 0 $((kptlines-1))); do
	wt=$(echo $kpthighsym | awk '{print $((row+2)*4)}' row=$row)
	if [[ $wt != "start" && $row != $kptlines ]]; then
		for i in $(seq 1 $kptpathden); do
			kxprev=$(echo $kpthighsym | awk '{print $((row*4)+1)}' row=$row)
			kxnext=$(echo $kpthighsym | awk '{print $(((row+1)*4)+1)}' row=$row)
			kyprev=$(echo $kpthighsym | awk '{print $((row*4)+2)}' row=$row)
			kynext=$(echo $kpthighsym | awk '{print $(((row+1)*4)+2)}' row=$row)
			kzprev=$(echo $kpthighsym | awk '{print $((row*4)+3)}' row=$row)
			kznext=$(echo $kpthighsym | awk '{print $(((row+1)*4)+3)}' row=$row)
			kptpath="${kptpath}\n$(python -c "print $kxprev+$i*($kxnext-$kxprev)/float($kptpathden)") $(python -c "print $kyprev+$i*($kynext-$kyprev)/float($kptpathden)") $(python -c "print $kzprev+$i*($kznext-$kzprev)/float($kptpathden)") $kptwt"
		done
	else
		kptpath="${kptpath}\n$(echo $kpthighsym | awk '{print $(((row+1)*4)+1), $(((row+1)*4)+2), $(((row+1)*4)+3)}' row=$row) $kptwt"
	fi
done

#Use the same $filename.save output from SCF for Berry-phase calculations
if [[ $loop != *"p"* ]]; then cp -r ./save_temp/$filename.save ./; fi


echo "Calculating the band structure."
cat > $filename.pwbands.in << EOF
&CONTROL
		 calculation = 'bands',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'restart',
		  wf_collect = .true.,
		   verbosity = 'high'
 /

 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
                         nat = $nat,
                        ntyp = $ntyp,
			nbnd = $nbnd, 
                     ecutwfc = $ewfc,
                 occupations = 'fixed' ,
                     degauss = 0.00 ,
 /
 &ELECTRONS
                    conv_thr = 1.0e-6,
                 mixing_beta = 0.3 ,
 /
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_vc

K_POINTS crystal_b
$(echo -e "$kptpath")

CELL_PARAMETERS alat
$a1_i $a2_i $a3_i
$b1_i $b2_i $b3_i
$c1_i $c2_i $c3_i
EOF
if [[ $loop != *"p"* && ( ! -a $pwbandsfile || -z $(grep 'PWSCF  ' $pwbandsfile) ) ]]; then $runprefix pw.x < $filename.pwbands.in > $pwbandsfile; fi

#Get the runtime
pwbandstime=$(grep "PWSCF" $pwbandsfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')


echo "Arranging the k-points of the band structure."
cat > $filename.bands.in << EOF
&BANDS
	 prefix = '$filename'
	filband = '$filename.bands.dat'
/
EOF
bandsfile="$filename.bands.out"
#if [[ $loop != *"p"* && ( ! -a $bandsfile || -z $(grep 'BANDS  ' $bandsfile) ) ]]; then $runprefix pw.x < $filename.bands.in > $bandsfile; fi
if [[ $loop != *"p"* && ( ! -a $bandsfile || -z $(grep 'BANDS  ' $bandsfile) ) ]]; then bands.x < $filename.bands.in > $bandsfile; fi

#Get the runtime
bandstime=$(grep "BANDS" $bandsfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')

echo -e "\nTotal Runtime\t SCF: $scftime WALL\t NSCF: $nscftime WALL\t DOS: $dostime WALL\t PDOS: $pdostime WALL\tPWBANDS: $pwbandstime WALL\tBANDS: $bandstime WALL" >> $output

#Rearrange data of bands.dat.gnu because a different format is presented in v6.0 compared to v5.3
if [[ $loop = *"6.0"* ]]; then
	bandtemp="$filename.bands.dat.temp"
	banddata="$filename.bands.dat.gnu"
	bandloc=$(( ($kpttot - 1)*$kptwt + 2 ))
	bandstart=1
	bandend=$bandloc
	if [[ -a $bandtemp ]]; then rm $bandtemp; fi
	
	for bandnum in $(seq 1 $nbnd); do
		if [[ $((bandnum%2)) = 1 ]]; then
			echo "$(cat $banddata)" | sed -n "${bandstart},${bandend}p" >> $bandtemp
		else
			bandmult=$((nbnd - (bandnum-1) - 1))
			echo "$(tac $banddata)" | sed -n "$((bandstart + bandmult*bandloc)),$((bandend + bandmult*bandloc))p" >> $bandtemp
		fi
		bandstart=$((bandstart+bandloc))
		bandend=$((bandend+bandloc))
	done
	rm $banddata; mv $bandtemp $banddata
fi

#Plotting the density of states
cat > $filename.dos.py << EOF
import numpy as np
import matplotlib.pyplot as plt

print "\nPlotting the density of states."

filename='$filename'
plotrange='$plotrange'
dosplot='$dos_plot'
fermi=$scffermi

Data = np.loadtxt('%s.dos' % (filename), str)
xData = np.array(Data[1:, 0], float)-fermi
yData1 = np.array(Data[1:, 1], float)
plt.figure(num=1, figsize=(10, 8))

xmin, xmax = plotrange.split("to")

plt.xlabel('Energy (eV)', size=24)
plt.ylabel('Density of States', size=24)
plt.plot(xData, yData1, 'r', label='DOS')
plt.legend(frameon=False, loc='best')
plt.xlim(float(xmin), float(xmax))

if dosplot != []:
	plt.savefig('%s.DOS.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()

#yData2 = np.array(Data[1:, 2], float)
#plt.figure(num=1, figsize=(10, 8))

#plt.xlabel('Energy (eV)', size=24)
#plt.ylabel('Density of States', size=24)
#plt.plot(xData, yData1, 'r', label='DOS')
#plt.plot(xData, yData2, 'g--', label='PDOS')
#plt.legend(frameon=False, loc='best')
#plt.xlim(float(xmin), float(xmax))

#if pdosplot != []:
#	plt.savefig('%s.PDOS.%s.pdf' % (filename, plotrange), format='pdf')
#plt.close()
EOF

#Plotting the partial density of states
cat > $filename.pdos.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import os, fnmatch

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

print "\nPlotting the partial density of states."

filename='$filename'
plotrange='$plotrange'
ranges = plotrange.split(",")
pdosplot='$pdos_plot'
fermi=$scffermi
lineplot = ["r", "g--", "b--", "c--", "y--", "m--"]

Data = np.loadtxt('%s.pdos.pdos_tot' % (filename), str)
xData = np.array(Data[1:, 0], float)-fermi
yData = np.array(Data[1:, 1], float)

#------------------ATOM PROJECTED--------------------#
pdos_all = find('%s.pdos.pdos_atm*' % (filename), './')
pdos_atoms = []
pdos_atoms_sum = []
pdos_orbitals = []
pdos_orbitals_sum = []

for pdos_name in pdos_all:
	pdos_file = np.loadtxt(pdos_name, str)
	pdos_Data = np.array(pdos_file[1:, 1], float)
	#pdos_sum += np.array(pdos_Data[1:, 1], float)
	
	atom = pdos_name.split("(")[1].split(")")[0]
	if atom not in pdos_atoms:
		pdos_atoms.append(atom)
		pdos_atoms_sum.append(pdos_Data)
	else:
		pdos_atoms_sum[pdos_atoms.index(atom)] += pdos_Data

for range in ranges:
	plt.figure(num=1, figsize=(10, 8))
	xmin0, xmax0 = range.split("to")
	xmin, xmax = abs(xData-float(xmin0)).argmin(), abs(xData-float(xmax0)).argmin()
	ymax = max(yData[xmin:xmax])

	plt.xlabel('Energy (eV)', size=24)
	plt.ylabel('Density of States', size=24)
	plt.plot(xData, yData, lineplot[0], alpha=0.75, label='Total')
	
	for j in xrange(len(pdos_atoms)):
		ymax = max(ymax, max(pdos_atoms_sum[j][xmin:xmax]))
		plt.plot(xData, pdos_atoms_sum[j], lineplot[j+1], label=pdos_atoms[j])
	plt.legend(frameon=False, loc='best')
	plt.xlim(xData[xmin], xData[xmax])
	plt.ylim(0.0, ymax*1.2)
	
	if pdosplot != []:
		plt.savefig('%s.PDOS-atom.%s.pdf' % (filename, range), format='pdf')
	plt.close()

#----------------------------------------------------#

#---------------ORBITAL PROJECTED--------------------#
for pdos_name in pdos_all:
	pdos_file = np.loadtxt(pdos_name, str)
	pdos_Data = np.array(pdos_file[1:, 1], float)
	#pdos_sum += np.array(pdos_Data[1:, 1], float)
	
	orbital = pdos_name.split("(")[2].split(")")[0]
	if orbital not in pdos_orbitals:
		pdos_orbitals.append(orbital)
		pdos_orbitals_sum.append(pdos_Data)
	else:
		pdos_orbitals_sum[pdos_orbitals.index(orbital)] += pdos_Data

for range in ranges:
	plt.figure(num=1, figsize=(10, 8))
	xmin0, xmax0 = range.split("to")
	xmin, xmax = abs(xData-float(xmin0)).argmin(), abs(xData-float(xmax0)).argmin()
	ymax = max(yData[xmin:xmax])
	
	plt.xlabel('Energy (eV)', size=24)
	plt.ylabel('Density of States', size=24)
	plt.plot(xData, yData, lineplot[0], alpha=0.75, label='Total')
	
	for j in xrange(len(pdos_orbitals)):
		ymax = max(ymax, max(pdos_atoms_sum[j][xmin:xmax]))
		plt.plot(xData, pdos_orbitals_sum[j], lineplot[j+1], label=pdos_orbitals[j])
	plt.legend(frameon=False, loc='best')
	plt.xlim(xData[xmin], xData[xmax])
	plt.ylim(0.0, ymax*1.2)
	
	if pdosplot != []:
		plt.savefig('%s.PDOS-orbital.%s.pdf' % (filename, range), format='pdf')
	plt.close()
#----------------------------------------------------#
EOF

#Plotting the band structure
cat > $filename.bands.py << EOF
import numpy as np
import matplotlib.pyplot as plt

print "\nPlotting band structure."

filename='$filename'
output='$output'
reflevel=$scffermi
rydberg=$RytoeV
bandsplot='$bands_plot'

Data = np.loadtxt('%s.bands.dat.gnu' % (filename), float)
xData = np.array(Data[:, 0], float)
#yData = (np.array(Data[:, 1], float)*rydberg - reflevel)
yData = (np.array(Data[:, 1], float) - reflevel)
plt.figure(num=1, figsize=(14, 10))

plt.plot(xData, yData)
plt.xlim(0, xData[$bandloc-2]-0.01)
plt.ylim(-4, 6)

if bandsplot != []:
	plt.savefig('%s.Bands.pdf' % (filename), format='pdf')
plt.close()
EOF

if [[ ! -z $dos_plot ]]; then python $filename.dos.py; fi
if [[ ! -z $pdos_plot ]]; then python $filename.pdos.py; fi
if [[ ! -z $bands_plot ]]; then python $filename.bands.py; fi

done #for loop_arr
done #for material_arr
