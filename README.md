# MR-CDFT: Multi-Reference Covariant Density Functional Theory
MR-CDFT is a computational code for nuclear structure calculations for quadrupole-octupole deformed nuclei based on multi-reference covariant density functional theory.

---

## Features
- CDFT calculation

  Supports axial calculations within the covariant density functional theory for both even-even nuclei and odd-mass nuclei 

- Projection calculation

  Support only for even-even nuclei. Including the following projection: 
    - Particle-number projection (PNP)
    - Angular-momentum projection (AMP)
    - Parity projection (PP)

- Parallelization
    - MPI
    - OpenMP support for Projection calculations
- Cross-platform: Linux, macOS, Windows
- GUI
    Provide a graphical user interface (GUI) for personal PC. 

## Requirements
- Fortran 90/95 compiler (tested with gfortran (10.2.1) and ifort(2021.3.0))
- MPI 
- Intel MKL Library
- Optional: OpenMP for shared-memory parallelization
- Optional: Python 3 (for Graphical User Interface)

## Installation
1. Clone the repository:

    ```bash
    git clone https://github.com/liviler/MR_CDFT_f90.git
    cd MR_CDFT_f90
    ```

2. Compile the code:

    You can compile the code using either gfortran or Intel Fortran (ifort), depending on the compiler available on your system.
    * Using CMake 
        ```bash
            cmake --preset mpi-gfortran
            cmake --build --preset mpi-gfortran
        ```
    * Using gfortran
        ```bash
        make gfortran
        ```
    * Using ifort
        ```bash
        make ifort
        ```
        Before running the make command, ensure that GNU Make and the selected Fortran compiler are properly installed and available in your environment.

    After successful compilation, the executable, the executable `MRCDFT` will be generated in the `bin/` directory.

3. Adding to Environment Variables :

    To run `MRCDFT` from any directory, you need to add the executable path to your system’s PATH environment variable.
    * Linux / macOS
        1) Open your shell configuration file (e.g. `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`).
        2) Add the following line (replace the path with the actual location of the bin/ directory):
            ```bash
            export PATH="/full/path/to/MR_CDFT_f90/bin:$PATH"
            ```
        3) Reload the configuration file:
            ``` bash
            source ~/.bashrc
            ```
        4) Verify the installation:
            ```bash
            which MRCDFT
            ```
    * Windows
        1) Open System Properties → Advanced system settings → Environment Variables.
        2) Under User variables or System variables, select Path and click Edit.
        3) Add the full path to the bin/ directory containing MRCDFT.exe.
        4) Open a new command prompt and test:
            ```cmd
                MRCDFT
            ```
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

### Examples

To run a test calculation for $^{22}{\text{Ne}}$: 
```bash
cd examples/22Ne
bash run.sh
```

### Output Description
TODO: Modify the output file name.