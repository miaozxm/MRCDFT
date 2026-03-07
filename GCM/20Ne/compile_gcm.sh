#_______________________________________________________________________
#           compile the HWGCM code
#  projected density 
#_______________________________________________________________________
# load compiler to compile the code
# module load pgi
#_______________________________________________________________________
#  set path for the source file, workspace, etc.
#_______________________________________________________________________
export pathsource=..  #home/research/RPG/src/
export pathwork=./                    # directory where calculations are started 
export pathexec=${pathwork}/exec                       # directory where the exe files are stored

cd ${pathwork}
#_______________________________________________________________________
#  compile file 'makefile'
#_______________________________________________________________________
#FC = ifort -fpic -openmp -132 -O3 -funroll-loops
#FC = f90 -O3 -fpic -openmp -fastsse -Minform=inform   
#	rm -rf *.o
cat <<'EOF' > makefile
FC = gfortran  -fopenmp -mcmodel=large -ffixed-line-length-132 -ffree-line-length-none -fstack-arrays -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall # -Wl,--stack,1073741824
OBJ =  rasgcm.o

run: $(OBJ)
	${FC} -O2 -o run $(OBJ) 
	rm -rf *.o
clean:
	rm -rf *.out

EOF

#_______________________________________________________________________
#  parameter file 'dic.par' which specifies the parameters
#_______________________________________________________________________

cat <<'EOF' > dic.par 
c-----------------------------------------------------------------------
c     Parameter file for:
c-----------------------------------------------------------------------
      character(500) :: WFS_DIR 		
c---- maximal number for GFV
      parameter (   igfv  =     100 )
      parameter (  igfvbc =     20 )
c
      parameter (  jmax  =      15 )
      parameter (  jmax2 =  jmax+jmax  )
c
      parameter (   pi1  = 3.141592653589793d0 )
      parameter (   pi2  =       pi1+pi1 )
      parameter (   pi4  =       pi2+pi2 )
c---- max number of configurations
      parameter (    maxq =     30 )
c---- mesh points for (r, theta, phi)
      parameter (    ngr     =        16 )
!-------------------------------------------------------
EOF

#cp ${pathsource}/rasgcm4OA_${ver}.f ./rasgcm.f
cp ${pathsource}/rasgcm4EE.f ./rasgcm.f
#cp ${pathsource}/hwdens02.f ./hwdens.f
#mv makefile ${pathsource}/makefile
#mv dic.par ${pathsource}/dic.par
#cat hwdens.f >>hwgcm.f

if [[ -d exec ]] 
then
   echo 'exec already exists'
else
   mkdir exec
fi
#cd ${pathsource}
# begin to compile
echo 'beginning to compile'
echo '.......................' 
make 
echo '.......................' 
echo 'compilation is done'
mv run ${pathexec}/
rm makefile
rm dic.par
rm *.f
exit






