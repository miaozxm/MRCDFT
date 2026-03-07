#_______________________________________________________________________
#    This script is used to calculate the projected kernels
#    for GCM calculation with the mean-field wave functions
#_______________________________________________________________________
#!/bin/bash -l
#PBS -N RPG-PKC
#PBS -j oe
#PBS -o RPG-PKC.out
#PBS -e TPG-PKC.err
#PBS -l walltime=1:23:00:00

ELE=Ne #Mg
A=20 # 24
NUC=$A${ELE}   # the DIC code generates the files with nucleus name as element name
Nf=6
mq=6 #55 #11 # 15 # 5 #13 #13 #13 # 15 # 016 # 014 #21 #013 # 010       # number of mesh point in q-space
iam=0 # 1 # 0 # 0 #      # 0: 1D; 1: 3D
iphi=07 #07 #09 # 01 # 07 #09 #09 # 01  #     # number of meshpoints in gauge angle 
nbet=12 #16 # 14 #14 # 14 #14 # 14 # 10 # 14 # 14 # 14 # 01           # 01 -> pure PNP

iBlock=0 # 0 for even-even; 1 for odd-mass

 if [[ ${iam} -eq 0 ]] 
  then  
   fout=fort.1D.${iphi}.${nbet}
  else
   fout=fort.3D.${iphi}.${nbet}
 fi

#...................... path.sh
DATAPATH=../../../examples/20Ne/output
pathwork=../../GCM/20Ne                    # directory where calculations are started 

pathexec=${pathwork}/exec                      # directory where the exe files are stored

cd ${pathexec}
export GCM_FILES_DIR=$DATAPATH
touch $GCM_FILES_DIR/HFB.wfs
echo $GCM_FILES_DIR

 if [[ -f data ]]; then rm data ; fi
 if [[ -f betgam.dat ]]; then rm betgam.dat ; fi
#---------------------------------------------------------------------------------- 
# parameter sets that could be used in the beyond mean-field calculations
# Force    =  PC,PC-F1               ! Parameterset of the Lagrangian
# Force    =  PC,PC-PK1              ! Parameterset of the Lagrangian
# Force    =  PC,DD-PC1              ! Parameterset of the Lagrangian
#--------------------------------------------------------------------
cat <<EOF > data 
${ELE} $A                                                ! nucleus under consideration
isoa     =    0                                          ! 0(not odd-A); 1 (is odd-A)                 
n0f      =   ${Nf}                                          ! eMax: number of HO shells 
iswi(1)  =   ${iam}                                      ! AMP    : (0) 1DAMP (1) 3DAMP
iphi     =   ${iphi}                                     ! number of meshpoints in gauge angle 
nbet     =   ${nbet}                                     !  01 -> Pure PNP; others, PNP+AMP
Jmax     =    2                                          ! maximal spin value to be                 
kmax     =    4                                          ! number of states for each spin          
Zeta     = 5.00E-03 5.00E-03 5.00E-03 5.00E-03 5.00E-03  ! cutoff in norm            
EOF
####################################  
#  pay attention to the order

#   0.10  0.000
#   0.15  0.200
#   0.20  0.100
cat <<EOF > betgam.dat 
             ${mq}                   ! number of mesh-point in q-space  
    0.00  0.00   0.00
    0.10  0.00   0.00
    0.20  0.00   0.00
    0.30  0.00   0.00
    0.40  0.00   0.00
    0.50  0.00   0.00
EOF
cat <<EOF > betgam2.dat 
             ${mq}                   ! number of mesh-point in q-space  
    0.00  0.00   0.00
    0.10  0.00   0.00
    0.20  0.00   0.00
    0.30  0.00   0.00
    0.40  0.00   0.00
    0.50  0.00   0.00
EOF
cat <<EOF > bk_betgam.dat 
  0.18  0.00     0
  0.18  0.05     0
  0.18  0.10     0
  0.18  0.15     0
  0.18  0.20     0
  0.18  0.25     0
  0.18  0.30     0
  0.18  0.35     0
EOF
pwd
./run > ${fout}.out

cat ${fout}.out
echo hwgcm calculation is finished
echo ...done

