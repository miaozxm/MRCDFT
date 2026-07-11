## Usage
Run the program using:
```bash
MRCDFT -p para.dat -d b23.dat
```
`para.dat` contains the main input parameters, including the nuclear name and mass number, the number of oscillator shells, and other calculation settings. `b23.dat` contains the grid of quadrupole–octupole deformation constraints for the nuclei. And you can use `scripts/run.sh`  to generate these two files.

Note: If you encounter issues caused by memory limits during execution, you can temporarily remove the stack size limit on Linux by running the following command.
```bash
ulimit -s unlimited
```
For Windows systems, you can increase the stack size by specifying a larger stack allocation during compilation.