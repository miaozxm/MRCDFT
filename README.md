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
    git clone -b master https://github.com/liviler/MRCDFT.git
    ```
2.  [Install: Windows + ifx](doc/installation/windows_ifx.md) 


## Running a Test Calculation

After `MRCDFT.exe` has been successfully built, open a terminal and navigate to the root directory of the MRCDFT project.

#### Set the number of threads

If you are using a **Command Prompt (cmd)** terminal, run:

```cmd
set OMP_NUM_THREADS=4
set MKL_NUM_THREADS=4
```
These commands set the number of OpenMP and MKL threads used by each MPI process.

If you are using a **Power Shell** (PS) terminal, run:
```ps
$env:OMP_NUM_THREADS=4
$env:MKL_NUM_THREADS=4
```

If you are using a **Bash** shell, run
```bash
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
```
#### Run with multiple processes
To run a test calculation for $^{22}\mathrm{Ne}$, execute the test calculation with:
```bash
cd examples/22Ne
mpiexec -np 2 ../../bin/MRCDFT -p 22Ne_para.dat -d 22Ne_b23.dat
```
This example will run with 2 MPI processes, with each process using 4 threads.

See [USage](./doc/Usage.md) for more information.