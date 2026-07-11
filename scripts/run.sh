#!/bin/bash 
#PBS -N MRCDFT
#PBS -o MRCDFT.out
#PBS -e MRCDFT.err 
#PBS -l select=1:ncpus=24    

export OMP_NUM_THREADS=24

# ************************** setting path ***********************************************
MRCDFT_bin=../../bin
Work_path=./



# ************************** setting nucleus *********************************************
ELE=Ne 
A=22    # mass number
Nf=6    # No. of HO shells
NUC=${A}${ELE}



# ******* The following settings are only valid for odd-Z or odd-N nuclei ************************
# block type
iBlock=0  # 0: Non-blocking; 1: Block the given energy level; 2: Block according to K^\pi
# block level (valid when blockType=1)
BlockLeveln=00
BlockLevelp=00
# block Kpi (valid when iBlock=2)
Kn=0  # K of block neutron,0: no, 1: 1/2, 2: 3/2 ...  
Pin=1 # Pi of block neutron, 1: + , -1: -       
Kp=0  # K of block proton, 0: no,  1: 1/2, 2: 3/2 ...  
Pip=1 # Pi of block proton, 1: + , -1: -       
# block method
bMethod=2 #  1: blocking -> convergence;   2:convergence -> block; 3: convergence -> block -> convergence


# **********************  The following setting are set for projection ******************
Jmax=6
nphi=05       # number of meshpoints for gauge angles
nbeta=12
iss=1         # Symmetry in Kernels. 0: All ; 1: Triangular Matrix ; 2: Diagonal elements only
# ...... q1
q1s=1
q1e=-1
# ...... q2
q2s=1
q2e=-1


# **************************create working directory *******************************

run_path=$Work_path${NUC}_${Nf}
mkdir -p $run_path
cd $run_path
rm -rf ${run_path}/output/*


# ************************ generate the input file for MRCDFT calculation ************************
#--------------------------------------------------------------------
# parameter sets that could be used in the beyond mean-field calculations
# Force    =  PC,PC-F1               ! Parameterset of the Lagrangian
#  V0      =  308.000   321.000
# Force    =  PC,PC-PK1              ! Parameterset of the Lagrangian
#  V0      =  349.500   330.000      ! pairing strength for delta-force
#  V0      =  314.550   346.500      ! 
#--------------------------------------------------------------------
echo " create input file ...."
cat <<EOF > ${NUC}_para.dat 
n0f,n0b  = ${Nf}  8                    ! number of oscillator shells
b0       = -2.448                   ! oscillator parameter (fm) of basis (If it is less than zero, use the empirical formula.)
beta0    =  0.00                    ! deformation parameter of basis
betas    =  0.50                    ! deformation beta2 of W.S. potential
bet3s    =  0.00                    ! deformation beta3 of W.S. potential
maxi     =  400                     ! maximal number of iterations
xmix     =  0.50                    ! mixing parameter
inin     =  1                       ! 1 (calc. from beginning); 0 (read saved pot.) 
${ELE} $A                           ! nucleus under consideration
Ide      =  4  4                    ! Pairing control: 1. no  2. Frozen  3. G   4. Delta
Delta    =  0.000000  0.000000      ! Frozen Gaps (neutrons and protons)
Ga       =  0.000000  0.000000      ! Pairing-Constants GG = GA/A
Delta-0  =  2.000000  2.000000      ! Initial values for the Gaps
Vpair    =  308.000   321.000       ! pairing strength for delta force
Force    =  PC-F1
icstr    =  2                       ! Quadratic constraint (no 0; beta2 1; b2+b3 2)
cspr     =  10.00                   ! Spring constant
cmax     =  1.000                   ! cutoff for dE/db
iRHB     =  0                       ! 0: BCS; 1: RHB
iBlock   =  ${iBlock}               ! 0: Non-blocking; 1: Block the given energy level; 2: Block according to K^\pi
bln      =  ${BlockLeveln}          ! block level of Neutron (valid when iBlock=1)
blp      =  ${BlockLevelp}          ! block level of Proton  (valid when iBlock=1)
Kn       =  ${Kn}                   ! K of block neutron, 1: 1/2, 2: 3/2 ...  (valid when iBlock=2)
Pin      =  ${Pin}                  ! Pi of block neutron, 1: + , -1: -       (valid when iBlock=2)
Kp       =  ${Kp}                   ! K of block proton,  1: 1/2, 2: 3/2 ...  (valid when iBlock=2)
Pip      =  ${Pip}                  ! Pi of block proton, 1: + , -1: -        (valid when iBlock=2)
bMethod  =  ${bMethod}              ! block method, 1: blocking -> convergence;   2:convergence -> block; 3: convergence -> block -> convergence
Erot     =  1                       ! 1:Belyaev formula; 2: Nilsson formula; 3: Odd A formula
c-------------------------------------------------------------------
ProjType =  1                       ! 0 : only RMF 1: RMF+AMP 2: only AMP
AMP      =  1                       ! Angular Momentum Projection    : (0) no (1) 1DAMP (2) 3DAMP
PNP      =  1                       ! Particle Number Projection     : (0) no (1) yes
Kernels  =  $iss                    ! Symmetry in Kernels. 0: All ; 1: Triangular Matrix ; 2: Diagonal elements only
q1       =  $q1s   $q1e             ! Start and end indices for q1 parameters; end = -1 indicates that the end index is set to the total number of deformation parameters.
q2       =  $q2s   $q2e             ! Start and end indices for q2 parameters; end = -1 indicates that the end index is set to the total number of deformation parameters.
Jmax     =  $Jmax                   ! maximal spin value to be
icm      =  1                       ! cent-of-mass corr. 1: average  ; 2: HO approximation
nphi     =  ${nphi}                       ! number of meshpoints in gauge angles phi (odd numbers)
EulerSym =  1                       ! Symmetry of Euler angles (alpha, beta ,gamma ). 0: no, 1: Axially, 2: D2 
nalpha   =  1                       ! number of meshpoints in Euler angles alpha (even numbers) 
nbeta    =  ${nbeta}                      ! number of meshpoints in Euler angles beta (even numbers) 
ngamma   =  1                       ! number of meshpoints in Euler angles gamma (even numbers) 
lambda   =  3                       ! max lambda of 1B transition density
EOF

#-----------------------------------
# mesh points in deformation q-space
#----------------------------------- 
cat <<EOF > ${NUC}_b23.dat
    beta2 beta3  others
    0.20  0.30   0.00
    0.20  0.00   0.00
    0.10  0.30   0.00
    0.10  0.00   0.00
EOF

start_time=$(date +%s)
echo -e "\033[32m run ...\033[0m"

${MRCDFT_bin}/MRCDFT -p ${NUC}_para.dat -d ${NUC}_b23.dat

echo calculation is finished !
end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_minutes=$((execution_time / 60))
execution_time_seconds=$((execution_time % 60))
echo "Time cost : ${execution_time_minutes}min${execution_time_seconds}s"
echo Done!



