# Installation Guide (Windows + ifx)

This guide describes how to build the project on Windows using CMake, Visual Studio, and Intel oneAPI.

## Prerequisites

Before building the project, please install the following dependencies:

###  Required
1. **CMake**  
   Download and install CMake from:  

   https://cmake.org/download/

2. **Visual Studio**  
   Visual Studio is required to provide the Microsoft linker (link.exe) and Windows SDK, which ifx relies on for linking Fortran code into an executable.    If you have already installed Visual Studio, please verify that the version meets the requirements of the Intel oneAPI Toolkit.
   If not installed, download and install Visual Studio Community from:

   https://visualstudio.microsoft.com/zh-hans/downloads/

   During installation, make sure to select the following workload:
   - **Desktop development with C++**

   This workload provides the required MSVC compiler and build tools.

3. **Intel oneAPI Toolkit**  
   Intel oneAPI Toolkit provides the required Fortran compiler (`ifx`), Intel MPI Library, and Intel Math Kernel Library (MKL) needed for building the project. Download and install Intel oneAPI Toolkit from: 
    
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/oneapi-toolkit-download.html



>  [!Note]
> * Please remember the installation paths of Visual Studio and Intel oneAPI Toolkit, as they will be needed during the later compilation process.
> * After installing CMake, the environment variables are typically configured automatically. Open a terminal and run `cmake --version` to verify it works. If the command is not recognized, manually add CMake's bin directory (e.g., `C:\Program Files\CMake\bin`) to your system's PATH.
> * After installing the Intel oneAPI Toolkit, the environment variable is typically configured automatically. Open a terminal and run `mpiexec --version` to verify `mpiexec` works. If the command is not recognized, manually add the bin directory where mpiexec is located (e.g., `C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin`) to your system PATH.

### Optional
- NSIS (only for installer packaging)
- WiX Toolset (only for installer packaging)


## Environment Setup

After installation, open a **Command Prompt** (cmd) terminal and set the Visual Studio and Intel oneAPI environments:

```cmd
call "C:\Program Files (x86)\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```
> [!CAUTION]
> The installation paths may differ depending on your local setup. Please replace the paths above with the actual locations of your `vcvars64.bat` and `setvars.bat` files.

These two commands set the required compilation environment variables for the current terminal session.
- `vcvars64.bat` sets the Visual Studio C/C++ build environment, including the MSVC compiler, linker, Windows SDK, and related build tools.
- `setvars.bat` sets the Intel oneAPI environment, including the Intel Fortran compiler (`ifx`/`ifort`), Intel MPI Library, Intel MKL, and other Intel development tools.


## Build Instructions
Use the **Command Prompt** terminal in which the Visual Studio and Intel oneAPI environments have already been set, and navigate to the **root directory of the MRCDFT** project.

### Configure the Build Environment
First, verify that the compilation environment has been configured correctly by running:
```cmd
cmake --preset mpi-ifx
```
This command configures the project using the predefined CMake preset mpi-ifx. It detects the Intel Fortran compiler, MPI environment, MKL libraries, and generates the corresponding build files in the `build\mpi-ifx\` directory. 

If the `build\mpi-ifx\` directory already exists, it is recommended to delete it before running this command. Otherwise, cached CMake configuration files from previous builds may cause unexpected issues.

If an error occurs, it usually means that the current terminal environment is not configured correctly. In this case, please check whether Visual Studio and Intel oneAPI Toolkit have been installed properly, and make sure the environment initialization commands were executed successfully.

### Build the Project
If the configuration step completes successfully, the compilation environment is ready. You can then build the project by running:

```cmd
cmake --build --preset mpi-ifx
```
This command compiles the source code and generates the executable files according to the configuration defined in the mpi-ifx preset.

After the build completes successfully, the executable file `MRCDFT.exe` will be generated in the project's `bin\` directory.

> [!TIP]
> To run `MRCDFT.exe` from any directory in the command line, you need to add the directory containing `MRCDFT.exe` to your system `PATH` environment variable.

---
---

## Install and Package (Option)
After the project is successfully built, you can optionally install and package the software using CMake and CPack.

### Install to Local Directory
CMake provides an installation step that copies the compiled executable and required runtime files into a clean directory structure.
To install the project, run:
```bash
cmake --install build\mpi-ifx --prefix "C:\Program Files\MRCDFT"
```
> [!TIP]
> * The `--prefix` option specifies the installation destination.
> * If not specified, the default installation path will be `C:\Program Files\MRCDFT` (may require administrator privileges).

### Package the Project
To create a distributable package (ZIP archive or installer) of the installed project, use cpack (CMake’s packaging tool).
```bash
cpack --preset mpi-ifx
```
After execution, navigate to `build/mpi-ifx/packages`, You will find the following generated packages: 
* ZIP archive (portable distribution)
* NSIS installer (Windows installer)
* WIX installer

> [!NOTE]
> To avoid issues where the NSIS installer fails to modify the system PATH (e.g., warning: PATH too long, installer unable to modify PATH) due to long path limitations:
> Ensure that NSIS is already installed, then download the **large strings build**  from [Special_Builds](https://nsis.sourceforge.io/Special_Builds). Then replace the corresponding files in your existing NSIS installation directory with those from the downloaded zip.

## Workflow
You can execute the full configure-build–package pipeline in a single command:
```bash
cmake --workflow --preset mpi-ifx-to-package
```
This command automatically performs the complete workflow defined in the preset, including:

This single command will automatically execute the following steps in sequence:
* Configure – Run `cmake --preset mpi-ifx`
* Build – Run `cmake --build --preset mpi-ifx`
* Package – Run `cpack --preset mpi-ifx`

