#!/bin/bash
#SBATCH -p A

#Initial directory
cd BaNiO3/Piezo/6.0-smart/

#Run settings
path=$PATH
export OMP_NUM_THREADS=1
runprefix="mpirun -np 4 pw.x" #salloc -p A -n 4

#System settings and conversion factors
pseudodir="/home/sand/pang/Research/Pseudo/"
forc=9.72344907103d-5
ewfc=73.49861765 #Ry
kptden=2000 #k-point density per reciprocal atom
Avogadro=6.022140857e23 #atoms/mol
BohrtoAng=0.529177249 #Angstrom/bohr
Ang3tocm3=1.0e-24 #cm^3/Angstrom^3
Rybohr3toGPa=14710.5164168 #GPa/(Ry/bohr^3)
RytoeV=13.605698066 #eV/Ry
eVtoJ=1.60217733e-19 #J/eV
boltz=1.38064852e-23 #J/K

#Select the materials to study
material_arr=("InN") #"InN" "BaTiO3_o" "BaGeO3" "CaNiO3"  "BaCO3_o" "BaSiO3" "TlNiO3"
for material in ${material_arr[@]}; do

if [[ $material = 'BaNiO3' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ni   58.6934	Ni.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'LaOF' ]]; then
loop_arr=('6.0-smart-min-phase')
crys_sys="cubic"
point_grp="-43m"
atomic_species="
La   138.90547	La.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF 
F    18.998403	F.GGA-PBE-paw.UPF "

elif [[ $material = 'ZnO' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Zn   65.38	Zn.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'KNbO3_c' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="cubic"
point_grp="m-3m"
atomic_species="
K    39.0983	K.GGA-PBE-paw.UPF 
Nb   92.90638	Nb.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'KNbO3_t' ]]; then
loop_arr=('6.0-smart-min-phase')
crys_sys="tet_high"
point_grp="4mm"
atomic_species="
K    39.0983	K.GGA-PBE-paw.UPF 
Nb   92.90638	Nb.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'KNbO3_o' ]]; then
loop_arr=('6.0-smart-min-phase')
crys_sys="ortho"
point_grp="mm2"
atomic_species="
K    39.0983	K.GGA-PBE-paw.UPF 
Nb   92.90638	Nb.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'KNbO3_r' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="rhombo_high"
point_grp="3m"
atomic_species="
K    39.0983	K.GGA-PBE-paw.UPF 
Nb   92.90638	Nb.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'PbTiO3_c' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="cubic"
point_grp="m-3m"
atomic_species="
Pb   207.2	Pb.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'PbTiO3_t' ]]; then
loop_arr=('6.0-smart-min-phase')
crys_sys="tet_high"
point_grp="4mm"
atomic_species="
Pb   207.2	Pb.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'AlAgTe2' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="tet_high"
point_grp="-42m"
atomic_species="
Al   26.981539	Al.GGA-PBE-paw.UPF 
Ag   107.8682	Ag.GGA-PBE-paw.UPF
Te   127.6	Te.GGA-PBE-paw.UPF "

elif [[ $material = 'InAgTe2' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="tet_high"
point_grp="-42m"
atomic_species="
In   114.818	In.GGA-PBE-paw.UPF 
Ag   107.8682	Ag.GGA-PBE-paw.UPF
Te   127.6	Te.GGA-PBE-paw.UPF "

elif [[ $material = 'NaBrO3' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="cubic"
point_grp="23"
atomic_species="
Na   22.989769	Na.GGA-PBE-paw.UPF  
Br   79.904	Br.GGA-PBE-paw.UPF  
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'NiAsS' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="cubic"
point_grp="23"
atomic_species="
Ni   58.6934	Ni.GGA-PBE-paw.UPF  
As   74.9216	As.GGA-PBE-paw.UPF  
S    32.065	S.GGA-PBE-paw.UPF "

elif [[ $material = 'GaPt' ]]; then
loop_arr=('6.0-smart-min-phase')
crys_sys="cubic"
point_grp="23"
atomic_species="
Ga   69.723	Ga.GGA-PBE-paw.UPF  
Pt   195.084	Pt.GGA-PBE-paw.UPF "

elif [[ $material = 'GeRh' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="cubic"
point_grp="23"
atomic_species="
Ge   72.64	Ge.GGA-PBE-paw.UPF  
Rh   102.9055	Rh.GGA-PBE-paw.UPF "

elif [[ $material = 'WOF4' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="tet_low"
point_grp="4"
atomic_species="
W    183.84	W.GGA-PBE-paw.UPF 
O    15.999	O.GGA-PBE-paw.UPF
F    18.998403	F.GGA-PBE-paw.UPF "

elif [[ $material = 'NbNO' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase')
crys_sys="tet_high"
point_grp="4mm"
atomic_species="
Nb   92.90638	Nb.GGA-PBE-paw.UPF 
N    14.0067	N.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'VBO4' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="tet_low"
point_grp="-4"
atomic_species="
V    50.9415	V.GGA-PBE-paw.UPF 
B    10.811	B.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BAsO4' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="tet_low"
point_grp="-4"
atomic_species="
B    10.811	B.GGA-PBE-paw.UPF 
As   74.9216	As.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'AlAsO4' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="tet_low"
point_grp="-4"
atomic_species="
Al   26.981539	Al.GGA-PBE-paw.UPF 
As   74.9216	As.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'InPS4' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="tet_low"
point_grp="-4"
atomic_species="
In   114.818	In.GGA-PBE-paw.UPF 
P    30.973762	P.GGA-PBE-paw.UPF
S    32.065	S.GGA-PBE-paw.UPF "

elif [[ $material = 'BN' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
B    10.811	B.GGA-PBE-paw.UPF  
N    14.0067	N.GGA-PBE-paw.UPF "

elif [[ $material = 'AlN' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Al   26.981539	Al.GGA-PBE-paw.UPF  
N    14.0067	N.GGA-PBE-paw.UPF "

elif [[ $material = 'GaN' ]]; then
loop_arr=('6.0-smart-halfscan-phase' '6.0-smart-kgrid-halfscan-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Ga   69.723	Ga.GGA-PBE-paw.UPF  
N    14.0067	N.GGA-PBE-paw.UPF "

elif [[ $material = 'InN' ]]; then
loop_arr=('6.0-smart-firstscan-phase' '6.0-smart-lastscan-phase')
#loop_arr=('6.0-smart-kgrid-firstscan' '6.0-smart-kgrid-firstscan-phase' '6.0-smart-kgrid-lastscan' '6.0-smart-kgrid-lastscan-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
In   114.818	In.GGA-PBE-paw.UPF  
N    14.0067	N.GGA-PBE-paw.UPF "

elif [[ $material = 'TlN' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Tl   204.3833	Tl.GGA-PBE-paw.UPF  
N    14.0067	N.GGA-PBE-paw.UPF "

elif [[ $material = 'BaGeO3' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="ortho"
point_grp="222"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ge   72.64	Ge.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaTiO3_r' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="rhombo_high"
point_grp="3m"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaTiO3_o' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="ortho"
point_grp="mm2"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaTiO3_t' ]]; then
loop_arr=('6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="tet_high"
point_grp="4mm"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Ti   47.867	Ti.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaCO3_r' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="rhombo_high"
point_grp="3m"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
C    12.0107	C.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaCO3_o' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="ortho"
point_grp="mm2"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
C    12.0107	C.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'BaSiO3' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="ortho"
point_grp="222"
atomic_species="
Ba   137.327	Ba.GGA-PBE-paw.UPF 
Si   28.0855	Si.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'CaNiO3' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="hex"
point_grp="6mm"
atomic_species="
Ca   40.078	Ca.GGA-PBE-paw.UPF 
Ni   58.6934	Ni.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "

elif [[ $material = 'TlNiO3' ]]; then
loop_arr=('6.0-smart-min' '6.0-smart-min-phase' '6.0-smart-kgrid-min' '6.0-smart-kgrid-min-phase')
crys_sys="rhombo_high"
point_grp="3m"
atomic_species="
Tl   204.3833	Tl.GGA-PBE-paw.UPF  
Ni   58.6934	Ni.GGA-PBE-paw.UPF
O    15.999	O.GGA-PBE-paw.UPF "
fi

#Finished strain magnitudes for different directions

#BaNiO3
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#LaOF
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#KNbO3_c
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#KNbO3_t
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005* 0.00 0.005* 0.01*
#6: -0.01 -0.005 0.00 0.005 0.01
#*STOP 3

#KNbO3_o
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#2: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01* -0.005 0.00 0.005 0.01
#5: -0.005* 0.00 0.005 0.01*
#6: -0.01 -0.005 0.00 0.005 0.01
#*STOP 3

#KNbO3_r
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: 
#4: 
#6: 

#PbTiO3_c
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#PbTiO3_t
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01

#AlAgTe2
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: 0.00 0.005
#6: 

#InAgTe2
#6.0-smart
#1: 
#3: 
#4: 
#6: 

#NaBrO3
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#NiAsS
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#GaPt
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#GeRh
#6.0-smart, vcrelax*
#1: -0.01* -0.005* 0.00* 0.005* 0.01*
#4: -0.01* -0.005* 0.00* 0.005* 0.01*
#*STOP 3

#WOF4
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01

#NbNO
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01

#VBO4
#6.0-smart (vc-relax good)
#1: -0.005* 0.00 0.005* 0.01*
#3: 
#4: 
#6: 
#*STOP3

#BN
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#AlN
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-crys-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#GaN
#6.0-smart, vcrelax*
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid, vcrelax*
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#*STOP 3

#InN
#6.0-smart
#1: -0.01 -0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
#3: -0.01 -0.007* -0.006 -0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
#4: -0.01 -0.007 -0.006 -0.005 -0.004 -0.003* -0.002* -0.001** 0.00 0.001* 0.002* 0.003* 0.004 0.005 0.006 0.007 0.008 0.009 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: -0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009
#3: -0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009
#4: -0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001* 0.00 0.001* 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009
#*STOP2 scf
#**STOP2 rel scf

#TlN
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005* 0.00 0.005* 0.01
#6.0-smart-kgrid
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#*STOP2 scf

#BaGeO3
#6.0-smart
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 
#6.0-smart-kgrid
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 

#BaTiO3_r
#6.0-smart (primitive)
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01

#BaTiO3_o
#6.0-smart **STOP2 ***STOP3
#1: -0.01 -0.005 0.00 0.005 0.01
#2: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005*** 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01**
#5: -0.01 -0.005*** 0.00 0.005*** 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 

#BaTiO3_t
#6.0-smart
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#5: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-clamp
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005* 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6: -0.01 -0.005 0.00 0.005 0.01
#*STOP 3

#BaCO3_r
#6.0-smart, vcrelax*
#1: -0.01 -0.005 0.00 0.005 0.01
#3: -0.01 -0.005 0.00 0.005 0.01
#4: -0.01 -0.005 0.00 0.005 0.01
#6.0-smart-kgrid, vcrelax*
#1: -0.01 -0.005 0.00* 0.005 0.01
#3: -0.01 -0.005 0.00* 0.005 0.01
#4: -0.01 -0.005 0.00* 0.005 0.01
#*STOP 3

#BaCO3_o
#6.0-smart
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 
#6.0-smart-kgrid
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 

#BaSiO3
#6.0-smart
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 
#6.0-smart-kgrid
#1: 
#2: 
#3: 
#4: 
#5: 
#6: 

#CaNiO3
#6.0-smart
#1: 0.00 0.005
#3: 
#4: 
#6.0-smart-kgrid
#1: 
#3: 
#4: 

#TlNiO3
#6.0-smart
#1: 
#3: 
#4: 
#6: 
#6.0-smart-kgrid
#1: 
#3: 
#4: 
#6: 

#Save directory and material name
cd ../../../$material
shdir=$(pwd); cd Piezo/6.0-smart/ 
filename=$(echo $material | awk '{print tolower($1)}')

#Data gathering and plotting ranges
stress_dirs=(1 2 3 4 5 6)
polar_dirs=(1 2 3)
for loop in ${loop_arr[@]}; do

#Load modules for the specific version
if [[ $loop = *"5.4-smart-k812"* ]]; then
	cd ../5.4-smart-k812/
	export PATH=/home/sand/pang/espresso-5.4.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/1.10.0
	kpt=8
	kptplus=4
elif [[ $loop = *"5.4-smart"* ]]; then
	cd ../5.4-smart/
	export PATH=/home/sand/pang/espresso-5.4.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/1.10.0
	kpt=6
	kptplus=4
elif [[ $loop = *"5.4"* ]]; then
	cd ../5.4/
	export PATH=/home/sand/pang/espresso-5.4.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/1.10.0
	kpt=6
	kptplus=4
elif [[ $loop = *"6.0-smart-kgrid"* ]]; then
	cd ../6.0-smart-kgrid/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptplus=4
elif [[ $loop = *"6.0-smart-k812"* ]]; then
	cd ../6.0-smart-k812/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=8
	kptplus=4
elif [[ $loop = *"6.0-smart"* ]]; then
	cd ../6.0-smart/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=6
	kptplus=4
elif [[ $loop = *"6.0"* ]]; then
	cd ../6.0/
	export PATH=/home/sand/pang/qe-6.0/bin:$path
	module purge; module load gcc/4.9.2 mkl/2016 anaconda slurm openmpi/2.0.1
	kpt=6
	kptplus=4
fi

#Set the magnitude of the strains and plotting range
if [[ $loop = *"smart"* && $loop != *"p"* ]]; then
	if [[ $loop = *"min"* ]]; then
		strain_arr=(0.00 0.005 0.01 -0.005 -0.01)
		plotrange="-1to+1"
	elif [[ $loop = *"mid"* ]]; then
		strain_arr=(0.00 0.002 0.004 0.006 0.008 0.01 -0.002 -0.004 -0.006 -0.008 -0.01)
		plotrange="-1to+1"
	elif [[ $loop = *"halfscan"* ]]; then
		strain_arr=(0.00 0.001 0.002 0.003 0.004 0.005 -0.001 -0.002 -0.003 -0.004 -0.005)
		plotrange="-0.5to+0.5"
	elif [[ $loop = *"firstscan"* ]]; then
		strain_arr=(-0.002 -0.003 -0.004 -0.005 -0.006 -0.007 -0.008 -0.009 -0.01)
		plotrange="-1to-0.2"
	elif [[ $loop = *"lastscan"* ]]; then
		strain_arr=(0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01)
		plotrange="+0.2to+1"
	elif [[ $loop = *"scan"* ]]; then
		strain_arr=(0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 -0.001 -0.002 -0.003 -0.004 -0.005 -0.006 -0.007 -0.008 -0.009) #0.01 -0.01
		plotrange="-1to+1"
	fi
else
	if [[ $loop = *"min"* ]]; then
		strain_arr=(-0.01 -0.005 0.00 0.005 0.01)
		plotrange="-1to+1"
	elif [[ $loop = *"mid"* ]]; then
		strain_arr=(-0.01 -0.008 -0.006 -0.004 -0.002 0.00 0.002 0.004 0.006 0.008 0.01)
		plotrange="-1to+1"
	elif [[ $loop = *"halfscan"* ]]; then
		strain_arr=(-0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005)
		plotrange="-0.5to+0.5"
	elif [[ $loop = *"firstscan"* ]]; then
		strain_arr=(-0.01 -0.009 -0.008 -0.007 -0.006 -0.005 -0.004) # -0.003 -0.002
		plotrange="-1to-0.2"
	elif [[ $loop = *"lastscan"* ]]; then
		strain_arr=(0.004 0.005 0.006 0.007 0.008 0.009 0.01) # 0.002 0.003 
		plotrange="+0.2to+1"
	elif [[ $loop = *"scan"* ]]; then
		strain_arr=(-0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001 0.00 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009) #0.01 -0.01
		plotrange="-1to+1"
	fi
fi

#Set the directions to apply strain and which to plot
if [[ $loop = *"p-xx"* ]]; then
	strain_dirs=(1)
	elastdir_plot="[(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)]"
	piezodir_plot="[(1, 1), (2, 1), (3, 1)]"
	enerdir_plot=""
elif [[ $loop = *"p-yy"* ]]; then
	strain_dirs=(2)
	elastdir_plot="[(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (6, 2)]"
	piezodir_plot="[(1, 2), (2, 2), (3, 2)]"
	enerdir_plot=""
elif [[ $loop = *"p-zz"* ]]; then
	strain_dirs=(3)
	elastdir_plot="[(1, 3), (2, 3), (3, 3), (4, 3), (5, 3), (6, 3)]"
	piezodir_plot="[(1, 3), (2, 3), (3, 3)]"
	enerdir_plot=""
elif [[ $loop = *"p-yz"* ]]; then
	strain_dirs=(4)
	elastdir_plot="[(1, 4), (2, 4), (3, 4), (4, 4), (5, 4), (6, 4)]"
	piezodir_plot="[(1, 4), (2, 4), (3, 4)]"
	enerdir_plot=""
elif [[ $loop = *"p-zx"* ]]; then
	strain_dirs=(5)
	elastdir_plot="[(1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (6, 5)]"
	piezodir_plot="[(1, 5), (2, 5), (3, 5)]"
	enerdir_plot=""
elif [[ $loop = *"p-xy"* ]]; then
	strain_dirs=(6)
	elastdir_plot="[(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (6, 6)]"
	piezodir_plot="[(1, 6), (2, 6), (3, 6)]"
	enerdir_plot=""

else #if collecting data or plotting normally

	if [[ $crys_sys = 'tet_high' && $loop = *"w5"* ]]; then
		strain_dirs=(1 3 5 6)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (3, 3), (5, 5), (6, 6)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'cubic' ]]; then
		strain_dirs=(1 4)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (4, 4)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'hex' ]]; then
		strain_dirs=(1 3 4)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (3, 3), (4, 4)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'tet_high' ]]; then
		strain_dirs=(1 3 4 6)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (3, 3), (4, 4), (6, 6)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'tet_low' ]]; then
		strain_dirs=(1 3 4 6)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (6, 1), (3, 3), (4, 4), (6, 6)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'rhombo_high' ]]; then
		strain_dirs=(1 3 4)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (4, 1), (3, 3), (4, 4)]"; else elastdir_plot=""; fi
	elif [[ $crys_sys = 'rhombo_low' ]]; then
		strain_dirs=(1 3 4) #6, (6, 6)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (3, 3), (4, 4)]"; else elastdir_plot=""; fi
	else #'ortho'
		strain_dirs=(1 2 3 4 5 6)
		if [[ $loop = *"p"* ]]; then elastdir_plot="[(1, 1), (2, 1), (3, 1), (2, 2), (3, 2), (3, 3), (4, 4), (5, 5), (6, 6)]"; else elastdir_plot=""; fi
	fi

	if [[ $loop = *"p"* ]]; then
		enerdir_plot="[all]"
		if [[ $point_grp = '23' || $point_grp = '-43m' || $point_grp = '622' || $point_grp = '422' ]]; then
			piezodir_plot="[(1, 4)]"
		elif [[ $point_grp = '-6m2' ]]; then
			piezodir_plot="[(2, 1)]"
		elif [[ $point_grp = '-6' ]]; then
			piezodir_plot="[(1, 1), (2, 1)]"
		elif [[ $point_grp = '32' ]]; then
			piezodir_plot="[(1, 1), (1, 4)]"
		elif [[ $point_grp = '-42m' ]]; then
			piezodir_plot="[(1, 4), (3, 6)]"
		elif [[ $point_grp = '6mm' || $point_grp = '4mm' ]]; then
			piezodir_plot="[(3, 1), (3, 3), (2, 4)]"
		elif [[ $point_grp = '222' ]]; then
			piezodir_plot="[(1, 4), (2, 5), (3, 6)]"
		elif [[ $point_grp = '3m' ]]; then
			piezodir_plot="[(2, 1), (3, 1), (3, 3), (2, 4)]"
		elif [[ $point_grp = '6' || $point_grp = '4' ]]; then
			piezodir_plot="[(3, 1), (3, 3), (1, 4), (2, 4)]"
		elif [[ $point_grp = '-4' ]]; then
			piezodir_plot="[(3, 1), (1, 4), (2, 4), (3, 6)]"
		elif [[ $point_grp = 'mm2' ]]; then
			piezodir_plot="[(3, 1), (3, 2), (3, 3), (2, 4), (1, 5)]"
		else #'3'
			piezodir_plot="[(1, 1), (2, 1), (3, 1), (3, 3), (1, 4), (2, 4)]"
		fi
	else
		piezodir_plot=""
		enerdir_plot=""
	fi
fi

#Set output file
if [[ $loop = *"p"* ]]; then output="$filename.$loop.$plotrange.txt"; else output="$filename.$loop.txt"; fi
echo -n "" > $output


#Array initializations
stress_plot=()
energy_plot=()
polar_plot=()
phase_plot=()
phase_i_plot=()
phase_e_plot=()
strain_plot="[" #c_plot="["
stressdir_plot="["
straindir_plot="["
polardir_plot="["
strainlen=${#strain_dirs[@]}
stresslen=${#stress_dirs[@]}
epsilen=${#strain_arr[@]}
polarlen=${#polar_dirs[@]}
for i in $(seq 1 $strainlen); do
   energy_plot+=("[")
   for j in $(seq 1 $stresslen); do
      stress_plot+=("[")
   done
   for k in $(seq 1 $polarlen); do
      polar_plot+=("[")
      phase_plot+=("[")
      phase_i_plot+=("[")
      phase_e_plot+=("[")
   done
done

#Material parameters
rundir=$(pwd); cd $shdir
cif="Standard/*.cif"
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
if [[ $loop != *"p"* && ( ! -a $convfile || -z $(grep 'PWSCF  ' $convfile) ) ]]; then $runprefix pw.x < $filename.vcrelax.in > $convfile; fi

alat=$(grep -n "CELL_PARAMETERS" $convfile | tail -1 | cut -f1 -d:)
read alat a1_i a2_i a3_i b1_i b2_i b3_i c1_i c2_i c3_i <<<$(echo $(sed -n "$((alat)),$((alat+3))p" $convfile) | awk '{gsub(")",""); print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}')
a_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($a1_i**2) + abs($a2_i**2) + abs($a3_i**2))")
b_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($b1_i**2) + abs($b2_i**2) + abs($b3_i**2))")
c_i=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($c1_i**2) + abs($c2_i**2) + abs($c3_i**2))")
alpha_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($b1_i*$c1_i + $b2_i*$c2_i + $b3_i*$c3_i)/($b_i*$c_i)))")
beta_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($c1_i*$a1_i + $c2_i*$a2_i + $c3_i*$a3_i)/($c_i*$a_i)))")
gamma_i=$(python -c "from math import acos, degrees; print degrees(acos(($BohrtoAng*$alat)**2*($a1_i*$b1_i + $a2_i*$b2_i + $a3_i*$b3_i)/($a_i*$b_i)))")
pos_i=$(grep -n "ATOMIC_POSITIONS" $convfile | tail -1 | cut -f1 -d:)
pos_i=$(sed -n "$((pos_i+1)),$((pos_i+nat))p" $convfile)
mass_i=0.0
for row in $(seq 0 $(echo $atomic_species | awk '{print (NF/3 -1)}')); do
	mass_i=$(python -c "print $mass_i + $(echo $atomic_species | awk '{print $(row*3 + 2)}' row=$row)*$(grep -o "$(echo $atomic_species | awk '{print $(row*3 + 1)}' row=$row)" <<< "$pos_i" | wc -l)/$Avogadro")
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

#Go through each direction in strain_dirs
for dir in $(seq 0 $((strainlen-1))); do

#Reset plots
strain_plot="[" #c_plot="["


#Go through each magnitude in strain_arr
for delta in $(seq 0 $((epsilen-1))); do

#Apply deformation tensor
if [ ${strain_dirs[$dir]} = 1 ]; then
a1=$(python -c "print $a1_i*(1+${strain_arr[$delta]})")
a2=$a2_i
a3=$a3_i
b1=$(python -c "print $b1_i*(1+${strain_arr[$delta]})")
b2=$b2_i
b3=$b3_i
c1=$(python -c "print $c1_i*(1+${strain_arr[$delta]})")
c2=$c2_i
c3=$c3_i
elif [ ${strain_dirs[$dir]} = 2 ]; then
a1=$a1_i
a2=$(python -c "print $a2_i*(1+${strain_arr[$delta]})")
a3=$a3_i
b1=$b1_i
b2=$(python -c "print $b2_i*(1+${strain_arr[$delta]})")
b3=$b3_i
c1=$c1_i
c2=$(python -c "print $c2_i*(1+${strain_arr[$delta]})")
c3=$c3_i
elif [ ${strain_dirs[$dir]} = 3 ]; then
a1=$a1_i
a2=$a2_i
a3=$(python -c "print $a3_i*(1+${strain_arr[$delta]})")
b1=$b1_i
b2=$b2_i
b3=$(python -c "print $b3_i*(1+${strain_arr[$delta]})")
c1=$c1_i
c2=$c2_i
c3=$(python -c "print $c3_i*(1+${strain_arr[$delta]})")
elif [ ${strain_dirs[$dir]} = 4 ]; then
a1=$a1_i
a2=$a2_i
a3=$(python -c "print $a3_i+${strain_arr[$delta]}*$a2_i")
b1=$b1_i
b2=$b2_i
b3=$(python -c "print $b3_i+${strain_arr[$delta]}*$b2_i")
c1=$c1_i
c2=$c2_i
c3=$(python -c "print $c3_i+${strain_arr[$delta]}*$c2_i")
elif [ ${strain_dirs[$dir]} = 5 ]; then
a1=$a1_i
a2=$a2_i
a3=$(python -c "print $a3_i+${strain_arr[$delta]}*$a1_i")
b1=$b1_i
b2=$b2_i
b3=$(python -c "print $b3_i+${strain_arr[$delta]}*$b1_i")
c1=$c1_i
c2=$c2_i
c3=$(python -c "print $c3_i+${strain_arr[$delta]}*$c1_i")
elif [ ${strain_dirs[$dir]} = 6 ]; then
a1=$a1_i
a2=$(python -c "print $a2_i+${strain_arr[$delta]}*$a1_i")
a3=$a3_i
b1=$b1_i
b2=$(python -c "print $b2_i+${strain_arr[$delta]}*$b1_i")
b3=$b3_i
c1=$c1_i
c2=$(python -c "print $c2_i+${strain_arr[$delta]}*$c1_i")
c3=$c3_i
fi

echo -e "alat: $alat\na1: $a1\ta2: $a2\ta3: $a3\nb1: $b1\tb2: $b2\tb3: $b3\nc1: $c1\tc2: $c2\tc3: $c3\n" >> $output
#c=$(python -c "from math import sqrt; print $BohrtoAng*$alat*sqrt(abs($c1**2) + abs($c2**2) + abs($c3**2))"); c_plot+="$c,"


echo -e "\nComputing the piezoelectric coefficient with ${strain_arr[$delta]} strain along direction ${strain_dirs[$dir]}."

#--------------------------------START SMART STRAIN--------------------------------#
if [[ $loop = *"smart"* && $loop != *"p"* ]]; then

#Check if strain 0.00 has been finished, otherwise run it first
cat > ${filename}.relax.in << EOF
&CONTROL
		 calculation = 'relax',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'from_scratch',
	       forc_conv_thr = $forc,
		     tstress = .true.,
		     tprnfor = .true.
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
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
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_i

K_POINTS automatic 
  $kpta $kptb $kptc   0 0 0 

CELL_PARAMETERS alat
$a1_i $a2_i $a3_i
$b1_i $b2_i $b3_i
$c1_i $c2_i $c3_i
EOF
rel_0file="$filename.${strain_dirs[$dir]}.0.00.relax.out"
if [[ ! -a $rel_0file || -z $(grep 'PWSCF  ' $rel_0file) ]]; then 
	echo "Applying zero strain."
	echo -e "Strain 0.00 along direction ${strain_dirs[$dir]}\n" >> $output
	$runprefix pw.x < $filename.relax.in > $rel_0file
fi
strain_0=0.00
pos_0=$(grep -n "ATOMIC_POSITIONS" $rel_0file | tail -1 | cut -f1 -d:)
pos_0=$(sed -n "$((pos_0+1)),$((pos_0+nat))p" $rel_0file)

#Set pos_prev and strain_prev
if [[ -z ${pos_prev+x} ]]; then
	echo -e "\nInitialize with data from zero run."
	pos_prev=$pos_0
	strain_prev=$strain_0
else
	echo -e "\nInitialize with data from previous run."
	pos_prev=$pos
	strain_prev=${strain_arr[$delta-1]}
fi

#Decide initial parameters for relaxation
if [[ $(python -c "print abs(${strain_arr[$delta]} - $strain_prev) <= abs(${strain_arr[$delta]} - $strain_0)") = True ]]; then
	echo "Using data of previous strain run."
	pos_curr=$pos_prev
	strain_curr=$strain_prev
else
	echo "Using data of zero strain run."
	pos_curr=$pos_0
	strain_curr=$strain_0
fi

echo -e "Performing smart strain at strain=${strain_arr[$delta]} with atomic positions from strain=$strain_curr.\n"
echo -e "Performing smart strain at strain=${strain_arr[$delta]} with atomic positions from strain=$strain_curr.\n" >> $output

#For normal case withot smart strain
else 
	pos_curr=$pos_i

fi
#---------------------------------END SMART STRAIN---------------------------------#

echo "Applying strain."
echo -e "Strain ${strain_arr[$delta]} along direction ${strain_dirs[$dir]}\n" >> $output
strain_plot+="${strain_arr[$delta]},"

cat > ${filename}.relax.in << EOF
&CONTROL
		 calculation = 'relax',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
		restart_mode = 'from_scratch',
	       forc_conv_thr = $forc,
		     tstress = .true.,
		     tprnfor = .true.
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
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
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos_curr

K_POINTS automatic 
  $kpta $kptb $kptc   0 0 0 

CELL_PARAMETERS alat
$a1 $a2 $a3
$b1 $b2 $b3
$c1 $c2 $c3
EOF
relfile="$filename.${strain_dirs[$dir]}.${strain_arr[$delta]}.relax.out"
if [[ $loop != *"p"* && ( ! -a $relfile || -z $(grep 'PWSCF  ' $relfile) ) && $output != *"berry-only"* ]]; then $runprefix pw.x < $filename.relax.in > $relfile; fi

#Get the optimized lattice parameters
pos=$(grep -n "ATOMIC_POSITIONS" $relfile | tail -1 | cut -f1 -d:)
pos=$(sed -n "$((pos+1)),$((pos+nat))p" $relfile)

#Get the stresses, runtime, and total energy
#Voigt notation for stress: {1:xx, 2:yy, 3:zz, 4:yz, 5:zx, 6:xy}
echo "Stresses for the relaxed cell (in GPa):" >> $output
str=$(grep -n "total   stress" $relfile | tail -1 | cut -f1 -d:)
str=$(sed -n "$((str+1)),$((str+3))p" $relfile)
read sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 <<<$(echo $str | awk '{print $1*conv, $8*conv, $15*conv, $9*conv, $3*conv, $2*conv}' conv=$Rybohr3toGPa)
echo -e "sigma_1: $sigma1\t sigma_2: $sigma2\t sigma_3: $sigma3\t sigma_4: $sigma4\t sigma_5: $sigma5\t sigma_6: $sigma6\n" >> $output
reltime=$(grep "PWSCF" $relfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')
reletot=$(grep '!' $relfile | tail -1 | awk '{print $(NF-1)}')
echo -e "Total energy of relaxed cell: $reletot Ry\n" >> $output
energy_plot[$dir]+="$reletot,"

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
		     tstress = .true.,
		     tprnfor = .true.,
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
                         nat = $nat,
                        ntyp = $ntyp,
                     ecutwfc = $ewfc,
                 occupations = 'fixed' ,
                     degauss = 0.00 ,
 /
 &ELECTRONS
                    conv_thr = 1.0e-12,
                 mixing_beta = 0.3 ,
 /
ATOMIC_SPECIES
$atomic_species

ATOMIC_POSITIONS crystal 
$pos

K_POINTS automatic 
  $kpta $kptb $kptc   0 0 0 

CELL_PARAMETERS alat
$a1 $a2 $a3
$b1 $b2 $b3
$c1 $c2 $c3
EOF
scffile="$filename.${strain_dirs[$dir]}.${strain_arr[$delta]}.scf.out"
berrylastfile="$filename.${strain_dirs[$dir]}.${strain_arr[$delta]}.berry-${polar_dirs[-1]}.out"
if [[ $loop != *"p"* && ( ( ! -a $scffile || -z $(grep 'PWSCF  ' $scffile) ) || ( ! -a $berrylastfile || -z $(grep 'PWSCF  ' $berrylastfile) ) ) ]]; then $runprefix pw.x < $filename.scf.in > $scffile; fi

#Get the stresses and runtime
#Voigt notation for stress: {1:xx, 2:yy, 3:zz, 4:yz, 5:zx, 6:xy}
echo "Stresses from self-consistent field calculation (in GPa):" >> $output
str=$(grep -n "total   stress" $scffile | tail -1 | cut -f1 -d:)
str=$(sed -n "$((str+1)),$((str+3))p" $scffile)
read sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 <<<$(echo $str | awk '{print $1*conv, $8*conv, $15*conv, $9*conv, $3*conv, $2*conv}' conv=$Rybohr3toGPa)
echo -e "sigma_1: $sigma1\t sigma_2: $sigma2\t sigma_3: $sigma3\t sigma_4: $sigma4\t sigma_5: $sigma5\t sigma_6: $sigma6\n" >> $output
sigma_arr=($sigma1 $sigma2 $sigma3 $sigma4 $sigma5 $sigma6)
for i in $(seq 0 $((stresslen-1))); do
stress_plot[$dir*$stresslen+$i]+="${sigma_arr[${stress_dirs[$i]}-1]},"
done
scftime=$(grep "PWSCF" $scffile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')

#Save $filename.save from SCF
if [[ $loop != *"p"* && ! -d "save_temp" ]]; then mkdir save_temp; fi
if [[ $loop != *"p"* ]]; then mv $filename.save ./save_temp; fi


#Go through each direction in polar_dirs
for gdir in $(seq 0 $((polarlen-1))); do
if [ ${polar_dirs[$gdir]} = 1 ]; then
   kptbig=$((kpta+kptplus))
   kptberry="  $kptbig $kptb $kptc   0 0 0 "
elif  [ ${polar_dirs[$gdir]} = 2 ]; then
   kptbig=$((kptb+kptplus))
   kptberry="  $kpta $kptbig $kptc   0 0 0 "
elif  [ ${polar_dirs[$gdir]} = 3 ]; then
   kptbig=$((kptc+kptplus))
   kptberry="  $kpta $kptb $kptbig   0 0 0 "
fi

#Use the same $filename.save output from SCF for Berry-phase calculations
if [[ $loop != *"p"* ]]; then cp -r ./save_temp/$filename.save ./; fi

echo "Calculating Berry-phase polarization for the optimized cell for gdir ${polar_dirs[$gdir]}."
cat > ${filename}.berry.in << EOF
&CONTROL
		 calculation = 'nscf',
		  pseudo_dir = '$pseudodir',
                      prefix = '$filename',
	       forc_conv_thr = $forc,
		      lberry = .true.,
			gdir = ${polar_dirs[$gdir]},
		      nppstr = $kptbig,
 /
 &SYSTEM
                       ibrav = 0,
		   celldm(1) = $alat,
                         nat = $nat,
                        ntyp = $ntyp,
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
$pos

K_POINTS automatic 
$kptberry

CELL_PARAMETERS alat
$a1 $a2 $a3
$b1 $b2 $b3
$c1 $c2 $c3
EOF
berryfile="$filename.${strain_dirs[$dir]}.${strain_arr[$delta]}.berry-${polar_dirs[$gdir]}.out"
if [[ $loop != *"p"* && ( ! -a $berryfile || -z $(grep 'PWSCF  ' $berryfile) ) ]]; then $runprefix pw.x < $filename.berry.in > $berryfile; fi

#Get the polarization and Berry phase along direction gdir and runtime
polar=$(grep -e 'polarization direction' $berryfile)
partial=$(echo $polar | awk '{print $(4+2*mult)}' mult=${polar_dirs[$gdir]})
polar=$(grep -e 'C/m^2' $berryfile | awk '{print $3*mult}' mult=$partial | sed s'/)$//')
echo -e "Berry Phase Polarization for direction ${polar_dirs[$gdir]}: $polar C/m^2" >> $output
polar_plot[$dir*$polarlen+($gdir)]+="$polar,"

#Total Berry Phase
read phase mod <<<$(grep "TOTAL PHASE" $berryfile | awk '{print $3, $5}' | sed s'/)$//')
echo -e "Total Berry Phase for direction ${polar_dirs[$gdir]}: $phase (mod $mod)" >> $output
phase=$(python -c "from math import pi; print $phase*(2*pi/$mod)")
phase_plot[$dir*$polarlen+($gdir)]+="$phase,"

#Ionic contribution to the Berry Phase
read phase_i mod_i <<<$(grep "Ionic Phase" $berryfile | awk '{print $3, $5}' | sed s'/)$//')
echo -e "Ionic Phase for direction ${polar_dirs[$gdir]}: $phase_i (mod $mod_i)" >> $output
phase_i=$(python -c "from math import pi; print $phase_i*(2*pi/$mod_i)")
phase_i_plot[$dir*$polarlen+($gdir)]+="$phase_i,"

#Electronic contribution to the Berry Phase
read phase_e mod_e <<<$(grep "Electronic Phase" $berryfile | awk '{print $3, $5}' | sed s'/)$//')
echo -e "Electronic Phase for direction ${polar_dirs[$gdir]}: $phase_e (mod $mod_e)" >> $output
phase_e=$(python -c "from math import pi; print $phase_e*(2*pi/$mod_e)")
phase_e_plot[$dir*$polarlen+($gdir)]+="$phase_e,"
done #for polarization directions

echo -e "\nTotal Runtime\t Relax: $reltime WALL\t SCF: $scftime WALL" >> $output

for pdir in $(seq 0 $((polarlen-1))); do
berryfile="$filename.${strain_dirs[$dir]}.${strain_arr[$delta]}.berry-${polar_dirs[$pdir]}.out"
berrytime=$(grep "PWSCF" $berryfile | tail -1 | awk -v FS="(CPU|WALL)" '{print $2}')
echo -e "Berry ${polar_dirs[$pdir]}: $berrytime WALL" >> $output
done #for pdir

echo -e "\n\n" >> $output
done #for loop across delta


#Close the arrays
for i in $(seq 0 $((stresslen-1))); do
   stress_plot[$dir*$stresslen+$i]=${stress_plot[$dir*$stresslen+$i]::-1}"],"
done
stress_plot[$dir*$stresslen]="["${stress_plot[$dir*$stresslen]}
stress_plot[$dir*$stresslen+$i]=${stress_plot[$dir*$stresslen+$i]::-1}"],"
energy_plot[$dir]=${energy_plot[$dir]::-1}"],"
straindir_plot+="${strain_dirs[$dir]},"
for j in $(seq 0 $((polarlen-1))); do
   polar_plot[$dir*$polarlen+$j]=${polar_plot[$dir*$polarlen+$j]::-1}"],"
   phase_plot[$dir*$polarlen+$j]=${phase_plot[$dir*$polarlen+$j]::-1}"],"
   phase_i_plot[$dir*$polarlen+$j]=${phase_i_plot[$dir*$polarlen+$j]::-1}"],"
   phase_e_plot[$dir*$polarlen+$j]=${phase_e_plot[$dir*$polarlen+$j]::-1}"],"
done
polar_plot[$dir*$polarlen]="["${polar_plot[$dir*$polarlen]}
polar_plot[$dir*$polarlen+$j]=${polar_plot[$dir*$polarlen+$j]::-1}"],"
phase_plot[$dir*$polarlen]="["${phase_plot[$dir*$polarlen]}
phase_plot[$dir*$polarlen+$j]=${phase_plot[$dir*$polarlen+$j]::-1}"],"
phase_i_plot[$dir*$polarlen]="["${phase_i_plot[$dir*$polarlen]}
phase_i_plot[$dir*$polarlen+$j]=${phase_i_plot[$dir*$polarlen+$j]::-1}"],"
phase_e_plot[$dir*$polarlen]="["${phase_e_plot[$dir*$polarlen]}
phase_e_plot[$dir*$polarlen+$j]=${phase_e_plot[$dir*$polarlen+$j]::-1}"],"

done #for loop across directions
stress_plot[0]="["${stress_plot[0]}
stress_plot[$dir*$stresslen+$i]=${stress_plot[$dir*$stresslen+$i]::-1}"]"
energy_plot[0]="["${energy_plot[0]}
energy_plot[$dir]=${energy_plot[$dir]::-1}"]"
polar_plot[0]="["${polar_plot[0]}
polar_plot[$dir*$polarlen+$j]=${polar_plot[$dir*$polarlen+$j]::-1}"]"
phase_plot[0]="["${phase_plot[0]}
phase_plot[$dir*$polarlen+$j]=${phase_plot[$dir*$polarlen+$j]::-1}"]"
phase_i_plot[0]="["${phase_i_plot[0]}
phase_i_plot[$dir*$polarlen+$j]=${phase_i_plot[$dir*$polarlen+$j]::-1}"]"
phase_e_plot[0]="["${phase_e_plot[0]}
phase_e_plot[$dir*$polarlen+$j]=${phase_e_plot[$dir*$polarlen+$j]::-1}"]"

for k in $(seq 0 $((stresslen-1))); do
stressdir_plot+="${stress_dirs[$k]},"
done

for l in $(seq 0 $((polarlen-1))); do
polardir_plot+="${polar_dirs[$l]},"
done

strain_plot=${strain_plot::-1}"]" #c_plot=${c_plot::-1}"]"
stressdir_plot=${stressdir_plot::-1}"]"
straindir_plot=${straindir_plot::-1}"]"
polardir_plot=${polardir_plot::-1}"]"

############################################################################################
######### Plot the stress, polarization, and total energy with respect to strain ###########
############################################################################################

#Plotting stress across different directions
cat > $filename.stressvsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc

print "\nPlotting stress across strains applied at different directions."

filename='$filename'
loop='$loop'
plotrange='$plotrange'
output='$output'
stresslen=$stresslen
strainlen=$strainlen
epsilen=$epsilen
strainplot=$strain_plot
stressdir=$stressdir_plot
straindir=$straindir_plot
elastdirplot=$elastdir_plot
C = np.empty((6, 6), float)
crys_sys='$crys_sys'
stability = False

#Plot settings
width=2
marker=16
legendsize=24
if crys_sys == 'cubic':
	lineplot = ["b--", "b--", "r--"]
	markerplot = ["bD", "b^", "rv"]
	colorplot = ["blue", "blue", "red"]
elif crys_sys in ['hex', 'tet_high']:
	lineplot = ["b--", "b--", "b--", "r--", "g--", "y--"]
	markerplot = ["bD", "b^", "bv", "rs", "g>", "yo"]
	colorplot = ["blue", "blue", "blue", "red", "green", "yellow"]
elif crys_sys in ['tet_low', 'rhombo_high']:
	lineplot = ["b--", "b--", "b--", "b--", "r--", "g--", "y--"]
	markerplot = ["bD", "b^", "bv", "bs", "r>", "go", "y*"]
	colorplot = ["blue", "blue", "blue", "blue", "red", "green", "yellow"]
elif crys_sys == 'rhombo_low':
	lineplot = ["b--", "b--", "b--", "b--", "b--", "r--", "g--"]
	markerplot = ["bD", "b^", "bv", "bs", "b>", "ro", "g*"]
	colorplot = ["blue", "blue", "blue", "blue", "blue", "red", "green"]
else: #'ortho'
	lineplot = ["b--", "b--", "b--", "r--", "r--", "g--", "y--", "m--", "c--"]
	markerplot = ["bD", "b^", "bv", "rs", "r>", "go", "y*", "m<", "cP"]
	colorplot = ["blue", "blue", "blue", "red", "red", "green", "yellow", "magenta", "cyan"]
plotarr = [[[] for j in xrange(strainlen)] for i in xrange(stresslen)]
for k in xrange(len(elastdirplot)):
	plotarr[stressdir.index(elastdirplot[k][0])][straindir.index(elastdirplot[k][1])] = [lineplot[k], markerplot[k], colorplot[k]]

xData0 = np.array($strain_plot, float)
xData1 = 0.5*xData0**2 + xData0
xData2 = xData0
yData = -np.array(${stress_plot[@]}, float)

#Record the data to plot
f=open(output, 'a')
f.write("\nDeformations xData1: %s" % xData1)
f.write("\nDeformations xData2: %s\n\n" % xData2)
for i in xrange(stresslen):
	for j in xrange(strainlen):
		f.write("Stress %s - Strain %s: %s\n" % (stressdir[i], straindir[j], yData[j][i]))
f.close

#Create the figure background
xmin, xmax = min(min(xData1), min(xData2)), max(max(xData1), max(xData2))
ymin, ymax = yData[0][0][0], yData[0][0][0]
for i in xrange(stresslen):
	for j in xrange(strainlen):
		if (stressdir[i], straindir[j]) in elastdirplot:
			ymin, ymax = min(ymin, min(yData[j][i])), max(ymax, max(yData[j][i]))
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.figure(num=1, figsize=(14, 9))
plt.plot((xmin-xpad, xmax+xpad), (0, 0), color='0')

f=open(output, 'a')
for j in xrange(strainlen):
	if crys_sys: #takes care of C[j][i] = C[i][j] for the independent constants
		start = j
	else:
		start = 0
	for i in xrange(start, stresslen):
		
		#Switch between Green-Lagrange strains for axial and shear types
		if straindir[j] in [1, 2, 3]:
			xData = xData1
		else:
			xData = xData2
		
		#Perform linear regression
		m, b, r, pvalue, stderr = sc.linregress(xData, yData[j][i])
		x_left, x_right = xmin, xmax
		y_left, y_right = m*x_left+b, m*x_right+b
		
		eq = '$ y\ =\ %.10fx + %.10f$\n$ R^2\ =\ %.10f$' % (m, b, r**2)
		print "stress %s - strain %s: %s" % (stressdir[i], straindir[j], eq)
		f.write("\nstress %s - strain %s: %s" % (stressdir[i], straindir[j], eq))
				
		if (stressdir[i], straindir[j]) in elastdirplot:
			
			#Plot the trendline
			plt.plot([x_left, x_right], [y_left, y_right], plotarr[i][j][0], linewidth=width)
			
			#Plot the data
			plt.plot(xData, yData[j][i], plotarr[i][j][1], label='$ S\mathrm{_%s\ by\ }E\mathrm{_%s\ only}$'%(stressdir[i], straindir[j]), markersize=marker, markeredgecolor=plotarr[i][j][2], linewidth=5)
		
		#Record the elastic stiffness constant
		C[stressdir[i]-1][straindir[j]-1] = m
f.close



#================Calculate the temperature difference between the state of least stress and ground state================#
if "temp" in loop:
	energyData=np.array(${energy_plot[@]}, float)
	energyGS = energyData[np.argmin(energyData)/epsilen][np.argmin(energyData)%epsilen]
	T_Data=(energyData-energyGS)*$eVtoJ/$boltz #in K
	stressRef=min(np.ravel(abs(yData)))
	
	locGS, locRef = [], []
	printGS = "\nThe ground state experiences the following stress with respect to the state of least stress:"
	printRef = "\n\nThe state of least stress is found at the following equivalent temperature/s from the ground state:"
	f=open(output, 'a')
	for i in xrange(stresslen):
		for j in xrange(strainlen):
			for k in xrange(epsilen):
				if abs(yData[j][i][k]) == stressRef:
					locRef.append((j, k, i))
	
	f.write("\n\n")
	for i in xrange(strainlen):
		for j in xrange(epsilen):
			if energyData[i][j] == energyGS:
				locGS.append((i, j))
			f.write("Strain Direction %s - Magnitude %s:\t %s Ry \t=\t %s K\n" % (straindir[i], strainplot[j], energyData[i][j], T_Data[i][j]))
	
	print printGS
	f.write(printGS)
	for pair in locGS:
		stressDiff = "\nAt %s strain along direction %s:\t %s GPa" % (strainplot[pair[1]], straindir[pair[0]], [allstress[pair[1]] for allstress in yData[pair[0]]]-stressRef)
		print stressDiff
		f.write(stressDiff)
	
	print printRef
	f.write(printRef)
	for pair in locRef:
		tempDiff = "\nAt %s strain along direction %s:\t %s K" % (strainplot[pair[1]], straindir[pair[0]], T_Data[pair[0]][pair[1]])
		print tempDiff
		f.write(tempDiff)
	f.close
#====================================================End for temp diff====================================================#


#============================================Apply symmetry for elastic tensor============================================#
print "\nC (GPa) before applying symmetry.\n%s" % np.array_str(C, precision=0, suppress_small=True)

#Applies symmetry to the elastic stiffness tensor and tests stability based on Mouhat (2014).
if crys_sys in ['cubic', 'hex', 'tet_high', 'tet_low', 'rhombo_high', 'rhombo_low']:
	C[1][1] = C[0][0]
	C[4][4] = C[3][3]
	
	if crys_sys == 'cubic':
		C[2][2] = C[0][0]
		C[2][1], C[2][0] = C[1][0], C[1][0]
		C[5][5] = C[3][3]
		
		if (C[0][0] > abs(C[1][0])) and (C[0][0] + 2*C[1][0]) and (C[3][3] > 0): #Cubic
			stability = True
			
	else:
		C[2][1] = C[2][0]
		
		if crys_sys in ['hex', 'rhombo_high', 'rhombo_low']:
			C[5][5] = 0.5*(C[0][0] - C[1][0])
			
			if crys_sys != 'hex':
				C[3][1], C[5][4] = -C[3][0], C[3][0]
				
				if crys_sys == 'rhombo_low':
					C[4][1], C[5][3] = -C[4][0], -C[4][0]
					
					if (C[0][0] > abs(C[1][0])) and (C[3][3] > 0) and (C[2][0]**2 < 0.5*C[2][2]*(C[0][0] + C[1][0])) and (C[3][0]**2 + C[4][0]**2 < C[3][3]*C[5][5]): #Rhombohedral low
						stability = True
						
				elif (C[0][0] > abs(C[1][0])) and (C[3][3] > 0) and (C[2][0]**2 < 0.5*C[2][2]*(C[0][0] + C[1][0])) and (C[3][0]**2 < C[3][3]*C[5][5]): #Rhombohedral high
					stability = True
					
			elif (C[0][0] > abs(C[1][0])) and (2*C[2][0]**2 < C[2][2]*(C[0][0] + C[1][0])) and (C[3][3] > 0): #Hexagonal
				stability = True
				
		elif crys_sys == 'tet_low':
			C[5][1] = -C[5][0]
			
			if (C[0][0] > abs(C[1][0])) and (2*C[2][0]**2 < C[2][2]*(C[0][0] + C[1][0])) and (C[3][3] > 0) and (2*C[5][0]**2 < C[5][5]*(C[0][0] - C[2][0])): #Tetragonal low
				stability = True
				
		elif (C[0][0] > abs(C[1][0])) and (2*C[2][0]**2 < C[2][2]*(C[0][0] + C[1][0])) and (C[3][3] > 0): #Tetragonal high
			stability = True
			
elif crys_sys == 'ortho':
	if (C[0][0] > 0) and (C[0][0]*C[1][1] > C[1][0]**2) and (C[0][0]*C[1][1]*C[2][2] + 2*C[1][0]*C[2][0]*C[2][1] - C[0][0]*C[2][1]**2 - C[1][1]*C[2][0]**2 - C[2][2]*C[1][0]**2 > 0) and (C[3][3] > 0) and (C[4][4] > 0) and (C[5][5] > 0): #Orthorhombic
		stability = True

#To apply symmetry for C[j][i] = C[i][j] for all constants
if crys_sys:
	for i in xrange(6):
		for j in xrange(i, 6):
			if abs(C[j][i]) > 0: C[i][j] = C[j][i]

print "\nStability for this crystal is %s.\n" % stability
f=open(output, 'a')
f.write("\n\nStability for this crystal is %s." % stability)
f.close
#====================================================End for symmetry====================================================#

plt.legend(frameon=False, loc='best', fontsize=legendsize)
plt.xlabel('$ \mathrm{Green-Lagrange\ Strain\ }E$', size=24)
plt.ylabel('$ \mathrm{Stress\ }S\ (\mathrm{GPa})$', size=24)
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-1*ypad, ymax+1.25*ypad)

#Format axis labels and tick labels with LaTeX
ax = plt.gca()
xticks = ax.get_xticks()
xlabels = []
for i in xrange(len(xticks)):
	xlabels.append('$%s$' % xticks[i])
ax.set_xticklabels(xlabels)
yticks = ax.get_yticks()
ylabels = []
for i in xrange(len(yticks)):
	ylabels.append('$%s$' % yticks[i])
ax.set_yticklabels(ylabels)
ax.tick_params(labelsize=20, length=10)#, pad=10)

if elastdirplot != []:
	plt.savefig('%s.StressStrain.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()

#Get elastic properties
if abs(C[5][5]) > 0: #if set([1, 2, 3, 4, 5, 6]) <= set(straindir) == set(stressdir):
	print "Deriving elastic properties."
	s = np.linalg.inv(C)
	K_V = ((C[0][0] + C[1][1] + C[2][2]) + 2*(C[0][1] + C[1][2] + C[2][0]))/9.0
	K_R = 1.0/((s[0][0] + s[1][1] + s[2][2]) + 2*(s[0][1] + s[1][2] + s[2][0]))
	G_V = ((C[0][0] + C[1][1] + C[2][2]) - (C[0][1] + C[1][2] + C[2][0]) + 3*(C[3][3] + C[4][4] + C[5][5]))/15.0
	G_R = 15.0/(4*(s[0][0] + s[1][1] + s[2][2]) - 4*(s[0][1] + s[1][2] + s[2][0]) + 3*(s[3][3] + s[4][4] + s[5][5]))
	K_VRH = (K_V + K_R)/2.0
	G_VRH = (G_V + G_R)/2.0
	AU = 5*(G_V/G_R) + (K_V/K_R) - 6
	mu = (3*K_VRH - 2*G_VRH)/(6*K_VRH + 2*G_VRH)
	print "C (GPa)\n%s\n\ns (10^{-12} Pa^{-1})\n%s\n\nK_V (GPa): %s\nK_R (GPa): %s\nG_V (GPa): %s\nG_R (GPa): %s\nK_{VRH} (GPa): %s\nG_{VRH} (GPa): %s\nA^U: %s\nmu: %s" % (np.array_str(C, precision=0, suppress_small=True), np.array_str(s/1e-3, precision=1, suppress_small=True), K_V, K_R, G_V, G_R, K_VRH, G_VRH, AU, mu)
	f=open(output, 'a')
	f.write("\n\nC (GPa)\n%s\n\ns (10^{-12} Pa^{-1})\n%s\n\nK_V (GPa): %.f\nK_R (GPa): %.f\nG_V (GPa): %.f\nG_R (GPa): %.f\nK_{VRH} (GPa): %.f\nG_{VRH} (GPa): %.f\nA^U: %.2f\nmu: %.2f" % (np.array_str(C, precision=1, suppress_small=True), np.array_str(s/1e-3, precision=1, suppress_small=True), K_V, K_R, G_V, G_R, K_VRH, G_VRH, AU, mu))
	f.close

EOF

#Plotting polarization across different directions
cat > $filename.polarvsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc

print "\nPlotting polarization vs. strain."

filename='$filename'
plotrange='$plotrange'
output='$output'
polarlen=$polarlen
strainlen=$strainlen
polardir=$polardir_plot
straindir=$straindir_plot
piezodirplot=$piezodir_plot
e = np.zeros((3, 6), float) #C = np.empty((polarlen, strainlen), float)
point_grp='$point_grp'

#Plot settings
width=2
marker=16
legendsize=24
if point_grp in ['23', '-43m', '622', '422', '-6m2', '32', '-42m', '6mm', '4mm', '222', 'mm2']:
	lineplot = ["b--", "r--", "g--", "c--", "y--", "m--"]
	markerplot = ["bD", "rs", "go", "c^", "yv", "m+"]
	colorplot = ["blue", "red", "green", "cyan", "yellow", "magenta"]
elif point_grp in ['-6', '3m']:
	lineplot = ["b--", "b--", "r--", "g--"]
	markerplot = ["bD", "bs", "ro", "g^"]
	colorplot = ["blue", "blue", "red", "green"]
elif point_grp in ['6', '4']:
	lineplot = ["b--", "r--", "g--", "g--"]
	markerplot = ["bD", "rs", "go", "g^"]
	colorplot = ["blue", "red", "green", "green"]
elif point_grp == '-4':
	lineplot = ["b--", "r--", "r--", "g--"]
	markerplot = ["bD", "rs", "ro", "g^"]
	colorplot = ["blue", "red", "red", "green"]
else: #'3'
	lineplot = ["b--", "b--", "b--", "r--", "g--", "g--"]
	markerplot = ["bD", "bs", "bo", "r^", "gv", "g+"]
	colorplot = ["blue", "blue", "blue", "red", "green", "green"]
plotarr = [[[] for j in xrange(strainlen)] for i in xrange(polarlen)]
for k in xrange(len(piezodirplot)):
	plotarr[polardir.index(piezodirplot[k][0])][straindir.index(piezodirplot[k][1])] = [lineplot[k], markerplot[k], colorplot[k]]

xData0 = np.array($strain_plot, float)
xData = xData0
#xData = 0.5*xData0**2 + xData0
yData = np.array(${polar_plot[@]}, float) #yData[straindir][polardir]

f=open(output, 'a')
f.write("\n\nStrain: %s\n\n" % xData)
for i in xrange(polarlen):
	for j in xrange(strainlen):
		f.write("Polarization %s - Strain %s: %s\n" % (polardir[i], straindir[j], yData[j][i]))
f.close


xmin, xmax = min(xData), max(xData)
ymin, ymax = yData[0][0][0], yData[0][0][0]
for i in xrange(polarlen):
	for j in xrange(strainlen):
		if (polardir[i], straindir[j]) in piezodirplot:
			ymin, ymax = min(ymin, min(yData[j][i])), max(ymax, max(yData[j][i]))
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.figure(num=1, figsize=(14, 9))

f=open(output, 'a')
for i in xrange(polarlen):
	for j in xrange(strainlen):
		m, b, r, pvalue, stderr = sc.linregress(xData, yData[j][i])
		x_left, x_right = xmin, xmax
		y_left, y_right = m*x_left+b, m*x_right+b
		eq = '$ y\ =\ %.10fx + %.10f$\n$ R^2\ =\ %.10f$' % (m, b, r**2)
		print "polarization %s - strain %s: %s" % (polardir[i], straindir[j], eq)
		f.write("\npolarization %s - strain %s: %s" % (polardir[i], straindir[j], eq))
		if (polardir[i], straindir[j]) in piezodirplot:
			plt.plot([x_left, x_right], [y_left, y_right], plotarr[i][j][0], linewidth=width)#, label=eq)
			plt.plot(xData, yData[j][i], plotarr[i][j][1], label='$ P\mathrm{_%s\ by\ }E\mathrm{_%s\ only}$'%(polardir[i], straindir[j]), markersize=marker, markeredgecolor=plotarr[i][j][2], linewidth=5)
		e[polardir[i]-1][straindir[j]-1] = m
f.close

print "\ne (C/m^2) before applying symmetry.\n%s" % np.array_str(e, precision=1, suppress_small=True)

#Applying symmetry to the piezoelectric tensor

if point_grp in ['23', '-43m', '-42m', '-4', '-6m2', '-6', '6mm', '4mm', '6', '4', '422', '622', '3m', '32', '3']:
	if point_grp in ['23', '-43m', '-42m', '-4']:
		e[1][4] = e[0][3]
		
		if point_grp in ['23', '-43m']:
			e[2][5] = e[0][3]
			
		elif point_grp == "-4":
			e[0][4] = -e[1][3]
			e[2][1] = -e[2][0]
			
	else:
		if point_grp in ['-6m2', '-6', '3m', '3']:
			e[1][1], e[0][5] = -e[1][0], e[1][0]
			
		if point_grp in ['-6', '32', '3']:
			e[0][1], e[1][5] = -e[0][0], -e[0][0]
			
		if point_grp in ['6', '4', '622', '422', '32', '3']:
			e[1][4] = -e[0][3]
			
		if point_grp in ['6mm', '4mm', '6', '4', '3m', '3']:
			e[2][1] = e[2][0]
			
			if point_grp == '3':
				e[0][4] = -e[1][3]
				
			else:
				e[0][4] = e[1][3]
				
print "\ne (C/m^2) after applying symmetry.\n%s" % np.array_str(e, precision=1, suppress_small=True)
f=open(output, 'a')
f.write("\n\ne (C/m^2)\n%s\n\n" % np.array_str(e, precision=5, suppress_small=True))
f.close

plt.legend(frameon=False, loc='best', fontsize=legendsize)
plt.xlabel('$ \mathrm{Strain\ }\epsilon$', size=24)
plt.ylabel('$ \mathrm{Berry-Phase Polarization\ }P\ (\mathrm{C/m}^2)$', size=24)
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-1*ypad, ymax+1.25*ypad)

#Format axis labels and tick labels with LaTeX
ax = plt.gca()
xticks = ax.get_xticks()
xlabels = []
for i in xrange(len(xticks)):
	xlabels.append('$%s$' % xticks[i])
ax.set_xticklabels(xlabels)
yticks = ax.get_yticks()
ylabels = []
for i in xrange(len(yticks)):
	ylabels.append('$%s$' % yticks[i])
ax.set_yticklabels(ylabels)
ax.tick_params(labelsize=20, length=10)#, pad=10)

if piezodirplot != []:
	plt.savefig('%s.PolarStrain.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()
EOF

#Plotting the total Berry phase across different directions
cat > $filename.phasevsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc
import scipy.constants as consts

print "\nPlotting total Berry phase vs. strain and the proper piezoelectric tensor."

filename='$filename'
plotrange='$plotrange'
output='$output'
polarlen=$polarlen
strainlen=$strainlen
crysdirlen=3
polardir=$polardir_plot
straindir=$straindir_plot
piezodirplot=$piezodir_plot
R = np.array([[$a1_i, $a2_i, $a3_i], [$b1_i, $b2_i, $b3_i], [$c1_i, $c2_i, $c3_i]], float)*$BohrtoAng*$alat/(1e10)
#dphi_deps = np.zeros((3, 6), float)
c_tilde = np.zeros((3, 6), float) #c = np.empty((polarlen, strainlen), float)
point_grp='$point_grp'
volume=$volume_i

#Plot settings
width=2
marker=16
legendsize=24
if point_grp in ['23', '-43m', '622', '422', '-6m2', '32', '-42m', '6mm', '4mm', '222', 'mm2']:
	lineplot = ["b--", "r--", "g--", "c--", "y--", "m--"]
	markerplot = ["bD", "rs", "go", "c^", "yv", "m+"]
	colorplot = ["blue", "red", "green", "cyan", "yellow", "magenta"]
elif point_grp in ['-6', '3m']:
	lineplot = ["b--", "b--", "r--", "g--"]
	markerplot = ["bD", "bs", "ro", "g^"]
	colorplot = ["blue", "blue", "red", "green"]
elif point_grp in ['6', '4']:
	lineplot = ["b--", "r--", "g--", "g--"]
	markerplot = ["bD", "rs", "go", "g^"]
	colorplot = ["blue", "red", "green", "green"]
elif point_grp == '-4':
	lineplot = ["b--", "r--", "r--", "g--"]
	markerplot = ["bD", "rs", "ro", "g^"]
	colorplot = ["blue", "red", "red", "green"]
else: #'3'
	lineplot = ["b--", "b--", "b--", "r--", "g--", "g--"]
	markerplot = ["bD", "bs", "bo", "r^", "gv", "g+"]
	colorplot = ["blue", "blue", "blue", "red", "green", "green"]
plotarr = [[[] for j in xrange(strainlen)] for i in xrange(polarlen)]
for k in xrange(len(piezodirplot)):
	plotarr[polardir.index(piezodirplot[k][0])][straindir.index(piezodirplot[k][1])] = [lineplot[k], markerplot[k], colorplot[k]]

xData0 = np.array($strain_plot, float)
xData = xData0
#xData = 0.5*xData0**2 + xData0
yData = np.array(${phase_plot[@]}, float) #yData[straindir][polardir]

f=open(output, 'a')
f.write("\n\nStrain: %s\n\n" % xData)
for i in xrange(polarlen):
	for j in xrange(strainlen):
		f.write("Berry Phase %s - Strain %s: %s\n" % (polardir[i], straindir[j], yData[j][i]))
f.close


xmin, xmax = min(xData), max(xData)
ymin, ymax = yData[0][0][0], yData[0][0][0]
for i in xrange(polarlen):
	for j in xrange(strainlen):
		if (polardir[i], straindir[j]) in piezodirplot:
			ymin, ymax = min(ymin, min(yData[j][i])), max(ymax, max(yData[j][i]))
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.figure(num=1, figsize=(14, 9))

f=open(output, 'a')
for i in xrange(polarlen):
	for j in xrange(strainlen):	
		m, b, r, pvalue, stderr = sc.linregress(xData, yData[j][i])
		x_left, x_right = xmin, xmax
		y_left, y_right = m*x_left+b, m*x_right+b
		eq = '$ y\ =\ %.10fx + %.10f$\n$ R^2\ =\ %.10f$' % (m, b, r**2)
		print "phase %s - strain %s: %s" % (polardir[i], straindir[j], eq)
		f.write("\nphase %s - strain %s: %s" % (polardir[i], straindir[j], eq))
		if (polardir[i], straindir[j]) in piezodirplot:
			plt.plot([x_left, x_right], [y_left, y_right], plotarr[i][j][0], linewidth=width)#, label=eq)
			plt.plot(xData, yData[j][i], plotarr[i][j][1], label='$ \phi\mathrm{_%s\ by\ }E\mathrm{_%s\ only}$'%(polardir[i], straindir[j]), markersize=marker, markeredgecolor=plotarr[i][j][2], linewidth=5)
		for alpha in xrange(crysdirlen):
			print "dphi/deps: %s\tR[alpha][i]: %s" % (m, R[alpha][polardir[i]-1])
			c_tilde[polardir[i]-1][straindir[j]-1] += m*R[alpha][polardir[i]-1]
f.close

c_tilde *= consts.e/(2*consts.pi*(volume/(1e10)**3))

print "\nc_tilde (C/m^2) before applying symmetry.\n%s" % np.array_str(c_tilde, precision=1, suppress_small=True)

#Applying symmetry to the piezoelectric tensor

if point_grp in ['23', '-43m', '-42m', '-4', '-6m2', '-6', '6mm', '4mm', '6', '4', '422', '622', '3m', '32', '3']:
	if point_grp in ['23', '-43m', '-42m', '-4']:
		c_tilde[1][4] = c_tilde[0][3]
		
		if point_grp in ['23', '-43m']:
			c_tilde[2][5] = c_tilde[0][3]
			
		elif point_grp == "-4":
			c_tilde[0][4] = -c_tilde[1][3]
			c_tilde[2][1] = -c_tilde[2][0]
			
	else:
		if point_grp in ['-6m2', '-6', '3m', '3']:
			c_tilde[1][1], c_tilde[0][5] = -c_tilde[1][0], c_tilde[1][0]
			
		if point_grp in ['-6', '32', '3']:
			c_tilde[0][1], c_tilde[1][5] = -c_tilde[0][0], -c_tilde[0][0]
			
		if point_grp in ['6', '4', '622', '422', '32', '3']:
			c_tilde[1][4] = -c_tilde[0][3]
			
		if point_grp in ['6mm', '4mm', '6', '4', '3m', '3']:
			c_tilde[2][1] = c_tilde[2][0]
			
			if point_grp == '3':
				c_tilde[0][4] = -c_tilde[1][3]
				
			else:
				c_tilde[0][4] = c_tilde[1][3]

print "\nc_tilde (C/m^2) after applying symmetry.\n%s" % np.array_str(c_tilde, precision=1, suppress_small=True)
f=open(output, 'a')
f.write("\n\nc_tilde (C/m^2)\n%s\n\n" % np.array_str(c_tilde, precision=5, suppress_small=True))
f.close


plt.legend(frameon=False, loc='best', fontsize=legendsize)
plt.xlabel('$ \mathrm{Strain\ }\epsilon$', size=24)
plt.ylabel('$ \mathrm{Berry-Phase\ }\phi$', size=24)
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-1*ypad, ymax+1.25*ypad)

#Format axis labels and tick labels with LaTeX
ax = plt.gca()
xticks = ax.get_xticks()
xlabels = []
for i in xrange(len(xticks)):
	xlabels.append('$%s$' % xticks[i])
ax.set_xticklabels(xlabels)
yticks = ax.get_yticks()
ylabels = []
for i in xrange(len(yticks)):
	ylabels.append('$%s$' % yticks[i])
ax.set_yticklabels(ylabels)
ax.tick_params(labelsize=20, length=10)#, pad=10)

if piezodirplot != []:
	plt.savefig('%s.PhaseStrain.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()
EOF

#Plotting the ionic Berry phase across different directions
cat > $filename.phase_i_vsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc

print "\nPlotting ionic Berry phase vs. strain."

filename='$filename'
plotrange='$plotrange'
output='$output'
polarlen=$polarlen
strainlen=$strainlen
polardir=$polardir_plot
straindir=$straindir_plot
piezodirplot=$piezodir_plot
point_grp='$point_grp'

#Plot settings
width=2
marker=16
legendsize=24
if point_grp in ['23', '-43m', '622', '422', '-6m2', '32', '-42m', '6mm', '4mm', '222', 'mm2']:
	lineplot = ["b--", "r--", "g--", "c--", "y--", "m--"]
	markerplot = ["bD", "rs", "go", "c^", "yv", "m+"]
	colorplot = ["blue", "red", "green", "cyan", "yellow", "magenta"]
elif point_grp in ['-6', '3m']:
	lineplot = ["b--", "b--", "r--", "g--"]
	markerplot = ["bD", "bs", "ro", "g^"]
	colorplot = ["blue", "blue", "red", "green"]
elif point_grp in ['6', '4']:
	lineplot = ["b--", "r--", "g--", "g--"]
	markerplot = ["bD", "rs", "go", "g^"]
	colorplot = ["blue", "red", "green", "green"]
elif point_grp == '-4':
	lineplot = ["b--", "r--", "r--", "g--"]
	markerplot = ["bD", "rs", "ro", "g^"]
	colorplot = ["blue", "red", "red", "green"]
else: #'3'
	lineplot = ["b--", "b--", "b--", "r--", "g--", "g--"]
	markerplot = ["bD", "bs", "bo", "r^", "gv", "g+"]
	colorplot = ["blue", "blue", "blue", "red", "green", "green"]
plotarr = [[[] for j in xrange(strainlen)] for i in xrange(polarlen)]
for k in xrange(len(piezodirplot)):
	plotarr[polardir.index(piezodirplot[k][0])][straindir.index(piezodirplot[k][1])] = [lineplot[k], markerplot[k], colorplot[k]]

xData0 = np.array($strain_plot, float)
xData = xData0
#xData = 0.5*xData0**2 + xData0
yData = np.array(${phase_i_plot[@]}, float) #yData[straindir][polardir]

f=open(output, 'a')
f.write("\n\nStrain: %s\n\n" % xData)
for i in xrange(polarlen):
	for j in xrange(strainlen):
		f.write("Ionic Phase %s - Strain %s: %s\n" % (polardir[i], straindir[j], yData[j][i]))
f.close


xmin, xmax = min(xData), max(xData)
ymin, ymax = yData[0][0][0], yData[0][0][0]
for i in xrange(polarlen):
	for j in xrange(strainlen):
		if (polardir[i], straindir[j]) in piezodirplot:
			ymin, ymax = min(ymin, min(yData[j][i])), max(ymax, max(yData[j][i]))
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.figure(num=1, figsize=(14, 9))

f=open(output, 'a')
for i in xrange(polarlen):
	for j in xrange(strainlen):
		m, b, r, pvalue, stderr = sc.linregress(xData, yData[j][i])
		x_left, x_right = xmin, xmax
		y_left, y_right = m*x_left+b, m*x_right+b
		eq = '$ y\ =\ %.10fx + %.10f$\n$ R^2\ =\ %.10f$' % (m, b, r**2)
		print "phase_i %s - strain %s: %s" % (polardir[i], straindir[j], eq)
		f.write("\nphase_i %s - strain %s: %s" % (polardir[i], straindir[j], eq))
		if (polardir[i], straindir[j]) in piezodirplot:
			plt.plot([x_left, x_right], [y_left, y_right], plotarr[i][j][0], linewidth=width)#, label=eq)
			plt.plot(xData, yData[j][i], plotarr[i][j][1], label='$ \phi\mathrm{_%s\ by\ }E\mathrm{_%s\ only}$'%(polardir[i], straindir[j]), markersize=marker, markeredgecolor=plotarr[i][j][2], linewidth=5)
f.close

plt.legend(frameon=False, loc='best', fontsize=legendsize)
plt.xlabel('$ \mathrm{Strain\ }\epsilon$', size=24)
plt.ylabel('$ \mathrm{Ionic\ Berry-Phase\ }\phi$', size=24)
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-1*ypad, ymax+1.25*ypad)

#Format axis labels and tick labels with LaTeX
ax = plt.gca()
xticks = ax.get_xticks()
xlabels = []
for i in xrange(len(xticks)):
	xlabels.append('$%s$' % xticks[i])
ax.set_xticklabels(xlabels)
yticks = ax.get_yticks()
ylabels = []
for i in xrange(len(yticks)):
	ylabels.append('$%s$' % yticks[i])
ax.set_yticklabels(ylabels)
ax.tick_params(labelsize=20, length=10)#, pad=10)

if piezodirplot != []:
	plt.savefig('%s.PhaseStrainIonic.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()
EOF

#Plotting the electronic Berry phase across different directions
cat > $filename.phase_e_vsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc

print "\nPlotting electronic Berry phase vs. strain."

filename='$filename'
plotrange='$plotrange'
output='$output'
polarlen=$polarlen
strainlen=$strainlen
polardir=$polardir_plot
straindir=$straindir_plot
piezodirplot=$piezodir_plot
point_grp='$point_grp'

#Plot settings
width=2
marker=16
legendsize=24
if point_grp in ['23', '-43m', '622', '422', '-6m2', '32', '-42m', '6mm', '4mm', '222', 'mm2']:
	lineplot = ["b--", "r--", "g--", "c--", "y--", "m--"]
	markerplot = ["bD", "rs", "go", "c^", "yv", "m+"]
	colorplot = ["blue", "red", "green", "cyan", "yellow", "magenta"]
elif point_grp in ['-6', '3m']:
	lineplot = ["b--", "b--", "r--", "g--"]
	markerplot = ["bD", "bs", "ro", "g^"]
	colorplot = ["blue", "blue", "red", "green"]
elif point_grp in ['6', '4']:
	lineplot = ["b--", "r--", "g--", "g--"]
	markerplot = ["bD", "rs", "go", "g^"]
	colorplot = ["blue", "red", "green", "green"]
elif point_grp == '-4':
	lineplot = ["b--", "r--", "r--", "g--"]
	markerplot = ["bD", "rs", "ro", "g^"]
	colorplot = ["blue", "red", "red", "green"]
else: #'3'
	lineplot = ["b--", "b--", "b--", "r--", "g--", "g--"]
	markerplot = ["bD", "bs", "bo", "r^", "gv", "g+"]
	colorplot = ["blue", "blue", "blue", "red", "green", "green"]
plotarr = [[[] for j in xrange(strainlen)] for i in xrange(polarlen)]
for k in xrange(len(piezodirplot)):
	plotarr[polardir.index(piezodirplot[k][0])][straindir.index(piezodirplot[k][1])] = [lineplot[k], markerplot[k], colorplot[k]]

xData0 = np.array($strain_plot, float)
xData = xData0
#xData = 0.5*xData0**2 + xData0
yData = np.array(${phase_e_plot[@]}, float) #yData[straindir][polardir]

f=open(output, 'a')
f.write("\n\nStrain: %s\n\n" % xData)
for i in xrange(polarlen):
	for j in xrange(strainlen):
		f.write("Electronic Phase %s - Strain %s: %s\n" % (polardir[i], straindir[j], yData[j][i]))
f.close


xmin, xmax = min(xData), max(xData)
ymin, ymax = yData[0][0][0], yData[0][0][0]
for i in xrange(polarlen):
	for j in xrange(strainlen):
		if (polardir[i], straindir[j]) in piezodirplot:
			ymin, ymax = min(ymin, min(yData[j][i])), max(ymax, max(yData[j][i]))
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.figure(num=1, figsize=(14, 9))

f=open(output, 'a')
for i in xrange(polarlen):
	for j in xrange(strainlen):
		m, b, r, pvalue, stderr = sc.linregress(xData, yData[j][i])
		x_left, x_right = xmin, xmax
		y_left, y_right = m*x_left+b, m*x_right+b
		eq = '$ y\ =\ %.10fx + %.10f$\n$ R^2\ =\ %.10f$' % (m, b, r**2)
		print "phase_e %s - strain %s: %s" % (polardir[i], straindir[j], eq)
		f.write("\nphase_e %s - strain %s: %s" % (polardir[i], straindir[j], eq))
		if (polardir[i], straindir[j]) in piezodirplot:
			plt.plot([x_left, x_right], [y_left, y_right], plotarr[i][j][0], linewidth=width)#, label=eq)
			plt.plot(xData, yData[j][i], plotarr[i][j][1], label='$ \phi\mathrm{_%s\ by\ }E\mathrm{_%s\ only}$'%(polardir[i], straindir[j]), markersize=marker, markeredgecolor=plotarr[i][j][2], linewidth=5)
f.close

plt.legend(frameon=False, loc='best', fontsize=legendsize)
plt.xlabel('$ \mathrm{Strain\ }\epsilon$', size=24)
plt.ylabel('$ \mathrm{Electronic\ Berry-Phase\ }\phi$', size=24)
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-1*ypad, ymax+1.25*ypad)

#Format axis labels and tick labels with LaTeX
ax = plt.gca()
xticks = ax.get_xticks()
xlabels = []
for i in xrange(len(xticks)):
	xlabels.append('$%s$' % xticks[i])
ax.set_xticklabels(xlabels)
yticks = ax.get_yticks()
ylabels = []
for i in xrange(len(yticks)):
	ylabels.append('$%s$' % yticks[i])
ax.set_yticklabels(ylabels)
ax.tick_params(labelsize=20, length=10)#, pad=10)

if piezodirplot != []:
	plt.savefig('%s.PhaseStrainElec.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()
EOF

#Plotting total energy across different directions
cat > $filename.etotvsstrain.py << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc

print "\nPlotting total energy vs. strain"

filename='$filename'
plotrange='$plotrange'
output='$output'
strainlen=$strainlen
epsilen=$epsilen
strainplot=$strain_plot
straindir=$straindir_plot
enerdirplot=$enerdir_plot

xData = np.array(strainplot, float)
yData = np.array(${energy_plot[@]}, float)
markerplot = ["rD-", "g^-", "bv-", "cs-", "yo-", "m8-"]

#Compute the temperature equivalent of the strained structures with respect to the ground state
energyGS = min(np.ravel(yData))
T_Data=(yData-energyGS)*$eVtoJ/$boltz #in K
printGS = "\nThe ground state energy of %s Ry can be found by applying" % energyGS
for i in xrange(strainlen):
	for j in xrange(epsilen):
		if yData[i][j] == energyGS:
			printGS += "\n%s strain along direction %s," % (strainplot[j], straindir[i])
printGS = printGS[:-1] + ".\n"
print printGS

f=open(output, 'a')
f.write("\n\n")
for i in xrange(strainlen):
	for j in xrange(epsilen):
		f.write("Strain Direction %s - Magnitude %s:\t %s Ry \t=\t %s K\n" % (straindir[i], strainplot[j], yData[i][j], T_Data[i][j]))
f.write(printGS)
f.close

plt.figure(num=1, figsize=(10, 6))
xmin, xmax = min(xData), max(xData)
ymin_i, ymax_i = np.argmin(yData), np.argmax(yData)
ymin, ymax = yData[ymin_i/epsilen][ymin_i%epsilen], yData[ymax_i/epsilen][ymax_i%epsilen]
xpad = (xmax-xmin)/5.0
ypad = (ymax-ymin)/5.0
plt.xlim(xmin-xpad, xmax+xpad)
plt.ylim(ymin-ypad, ymax+ypad)

for i in xrange(strainlen):
	plt.plot(xData, yData[i], markerplot[i], label='$ \mathrm{Strain\ %s}$'%(straindir[i]), markersize=8)
plt.legend(frameon=False, loc='best', fontsize='large')
plt.xlabel('$ \mathrm{Strain}$', size=18)
plt.ylabel('$ \mathrm{Total\ Energy}\ (\mathrm{Ry})$', size=18)

if enerdirplot != []:
	plt.savefig('%s.EnergyStrain.%s.pdf' % (filename, plotrange), format='pdf')
plt.close()
EOF

if [[ ! -z $elastdir_plot ]]; then python $filename.stressvsstrain.py; fi
if [[ ! -z $piezodir_plot ]]; then python $filename.polarvsstrain.py; fi
if [[ $loop = *"phase"* ]]; then python $filename.phasevsstrain.py; python $filename.phase_i_vsstrain.py; python $filename.phase_e_vsstrain.py; fi
if [[ ! -z $enerdir_plot ]]; then python $filename.etotvsstrain.py; fi

done #for loop_arr
done #for material_arr
