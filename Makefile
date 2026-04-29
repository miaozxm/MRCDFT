EXE_NAME = MRCDFT
ROOT_DIR := $(CURDIR)
EXE_DIR = ${ROOT_DIR}/bin
SRC_DIR = ${ROOT_DIR}/src
MOD_DIR = ${ROOT_DIR}/src/mod
OBJ_DIR = ${ROOT_DIR}/src/obj
SRC_FILE_PREFIX_1 = CDFT
SRC_FILE_PREFIX_2 = Proj
SRC_FILE_PREFIX_3 = GCM

SOURCES = $(wildcard $(SRC_DIR)/*.f90)
OBJECTS = $(patsubst $(SRC_DIR)/%.f90, ${OBJ_DIR}/%.o, $(SOURCES))

# export MKL_THREADING_LAYER = GNU

default:  gfortran #ifort

# compiled by gfortran
gfortran: FC = gfortran
gfortran: FFLAGS = -O3 -J ${MOD_DIR} -fopenmp -ffree-line-length-none
# gfortran: FFLAGS = -O3 -g -J ${MOD_DIR} -fopenmp -ffree-line-length-none
# gfortran: FFLAGS += -lmkl_rt
gfortran: FFLAGS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
# gfortran: FFLAGS = -O2 -J ${MOD_DIR}  -fopenmp -ffree-line-length-none  -fstack-arrays -Wl,--stack,1073741824 # use this line if you need to increase stack size
gfortran: printConfiguration ${EXE_NAME} printEndInformation

# compiled by gfortran with bounds checking and debug info
debug: FC = gfortran
# debug: FFLAGS = -std=legacy -g -fcheck=all -Wall -ffree-line-length-none -fopenmp -J ${MOD_DIR} 
# debug: FFLAGS = -std=legacy -g -fcheck=bounds -Wall -ffree-line-length-none -fopenmp -J ${MOD_DIR} 
debug: FFLAGS = -std=legacy -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -ffree-line-length-none -fopenmp -J ${MOD_DIR} 
debug: FFLAGS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
debug: printConfiguration ${EXE_NAME} printEndInformation

perf : FC = gfortran
perf : FFLAGS = -g -pthread -O3 -J ${MOD_DIR}  -fopenmp -ffree-line-length-none # -O3 -march=native -flto
perf : FFLAGS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
perf : printConfiguration ${EXE_NAME} printEndInformation

mpif90: FC = mpif90
mpif90: FFLAGS = -O3 -J ${MOD_DIR} -fopenmp -ffree-line-length-none
mpif90: FFLAGS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
mpif90: printConfiguration ${EXE_NAME} printEndInformation

# compiled by ifort
ifort: FC = ifort
ifort: FFLAGS = -O2 -module ${MOD_DIR} -qopenmp
ifort: FFLAGS += -mkl
ifort: printConfiguration ${EXE_NAME} printEndInformation

# compiled by ifort with bounds checking and debug info
debugifort: FC = ifort
debugifort: FFLAGS = -stand f95 -g -check bounds -warn all -traceback -fpp -qopenmp -module ${MOD_DIR}
debugifort: printConfiguration ${EXE_NAME} printEndInformation



printConfiguration:
	@echo "=============Compiling with ${FC}====================="
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "------------------------------------------------------"
	
printEndInformation:
ifeq ($(OS),Windows_NT)
	@echo -e "\033[32mCompilation Finished! \033[0m"
	@echo -e "MRCDFT path: \033[32m ${EXE_DIR}\\${EXE_NAME}.exe \033[0m"
	@echo ""
	@echo "To run '$(EXE_NAME)' from anywhere, add:"
	@echo -e "\033[33m ${EXE_DIR} \033[0m"
	@echo "to your PATH Environment Variables. Then start a new Command Prompt or PowerShell session."
else
	@echo -e "\033[32mCompilation Finished! \033[0m"
	@echo -e "MRCDFT path: \033[32m ${EXE_DIR}/${EXE_NAME} \033[0m"
	@echo ""
	@echo "To run '$(EXE_NAME)' from anywhere, add:"
	@echo -e "\033[33m  export PATH=\$$PATH:${EXE_DIR} \033[0m"
	@echo "to your ~/.bashrc or ~/.zshrc and reload the shell."
endif
${EXE_NAME}:${OBJECTS} | ${EXE_DIR}
	@echo "compiling $@ ......"
	${FC} ${FFLAGS} -o ${EXE_DIR}/${EXE_NAME} $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | ${OBJ_DIR}  ${MOD_DIR}
	@echo "compiling  $@ ......"
	$(FC) $(FFLAGS) -c $< -o $@ 

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

${MOD_DIR}:
	mkdir -p ${MOD_DIR}

${EXE_DIR}:
	mkdir -p ${EXE_DIR}

# Dependencies
${OBJ_DIR}/main.o : ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_main.o  ${OBJ_DIR}/${SRC_FILE_PREFIX_3}_main.o

${OBJ_DIR}/Globals.o : $(OBJ_DIR)/Constants.o ${OBJ_DIR}/Tools.o 
${OBJ_DIR}/Mathmethods.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

## CDFT Dependencies
${OBJ_DIR}/${SRC_FILE_PREFIX_1}_main.o : $(filter-out ${OBJ_DIR}/main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_main.o,  ${OBJ_DIR}/${SRC_FILE_PREFIX_3}_main.o ${OBJECTS})

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_field.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_inout.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_dirac_equation.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_force.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
						   			     ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_field.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_nucleus.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_basis.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_constraint.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_dirac_BCS.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o \
											  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_dirac_equation.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_dirac_equation.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_density.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_expectation.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_broyden.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o \
											${OBJ_DIR}/${SRC_FILE_PREFIX_1}_dirac_equation.o \
											${OBJ_DIR}/${SRC_FILE_PREFIX_1}_RHB_delta_field.o \
											${OBJ_DIR}/${SRC_FILE_PREFIX_1}_field.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_RHB_delta_field.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_RHB_equation.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o \
												 ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_broyden.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_expectation_rotation.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_1}_inout.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
					 			  		  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_expectation.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_expectation_rotation.o

## Proj Dependencies
${OBJ_DIR}/${SRC_FILE_PREFIX_2}_main.o : $(filter-out ${OBJ_DIR}/main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_main.o,  ${OBJ_DIR}/${SRC_FILE_PREFIX_3}_main.o ${OBJECTS})

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_kernel.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_inout.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_basis.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_mixed.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_Jsquare_Nsquare.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_electromagnetic_multipole.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_Energy.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_density.o \
										   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_eccentricity.o 

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_mixed.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o \
						                  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_basis.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_Energy.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_Jsquare_Nsquare.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_electromagnetic_multipole.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_density.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
											${OBJ_DIR}/${SRC_FILE_PREFIX_1}_basis.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_transition_density.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
									  				   ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_basis.o  \
													   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_electromagnetic_multipole.o
${OBJ_DIR}/${SRC_FILE_PREFIX_2}_eccentricity.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
													   ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_electromagnetic_multipole.o

${OBJ_DIR}/${SRC_FILE_PREFIX_2}_inout.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o\
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_inout.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_electromagnetic_multipole.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_eccentricity.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_transition_density.o 

## GCM Dependencies
${OBJ_DIR}/${SRC_FILE_PREFIX_3}_main.o : $(filter-out ${OBJ_DIR}/main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_1}_main.o ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_main.o,  ${OBJ_DIR}/${SRC_FILE_PREFIX_3}_main.o ${OBJECTS})

${OBJ_DIR}/${SRC_FILE_PREFIX_3}_HWG.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o ${OBJ_DIR}/Mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX_3}_observables.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX_3}_inout.o : ${OBJ_DIR}/Constants.o ${OBJ_DIR}/Globals.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX_2}_inout.o 


path:
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "SOURCES : ${SOURCES}"
	@echo "OBJECTS : ${OBJECTS}"

clean:
	rm -f ${EXE_DIR}/${EXE_NAME} $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

deepclean:
	rm -rf ${EXE_DIR} $(OBJ_DIR) $(MOD_DIR) ${SRC_DIR}/*.mod
