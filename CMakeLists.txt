cmake_minimum_required(VERSION 3.20)

project(MRCDFT LANGUAGES Fortran)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#  compile with Fortran 2008 standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
# set(CMAKE_Fortran_EXTENSIONS OFF)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0)
        message(FATAL_ERROR "GNU Fortran >= 9.0 is required.")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19.0)
        message(FATAL_ERROR "Intel Fortran >= 19.0 is required.")
    endif()
else()
    message(WARNING "Untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Create the executable target
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(EXE_DIR "${CMAKE_SOURCE_DIR}/bin")
set(MOD_DIR "${CMAKE_BINARY_DIR}/mod")
file(GLOB SRC_FILES CONFIGURE_DEPENDS
    "${SRC_DIR}/*.f90"
)
add_executable(MRCDFT ${SRC_FILES})
set_target_properties(MRCDFT PROPERTIES
    Fortran_MODULE_DIRECTORY "${MOD_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${EXE_DIR}"
)

# Set compiler options for different build types
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(MRCDFT PRIVATE
        -ffree-line-length-none
        $<$<CONFIG:Debug>:-g;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-Wall>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    target_compile_options(MRCDFT PRIVATE
        -free
        $<$<CONFIG:Debug>:-g;-check;all;-traceback;-warn;all>
    )
else()
    target_compile_options(MRCDFT PRIVATE
        $<$<CONFIG:Debug>:-g>
    )
endif()



# Find MPI, OpenMP and MKL libraries
find_package(MPI REQUIRED COMPONENTS Fortran)

find_package(OpenMP REQUIRED COMPONENTS Fortran)


set(MKL_LINK dynamic)
set(MKL_INTERFACE lp64)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MKL_THREADING gnu_thread)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(MKL_THREADING intel_thread)
else()
    set(MKL_THREADING sequential)
endif()
find_package(MKL CONFIG REQUIRED)

# Link against OpenMP and MKL libraries
target_link_libraries(MRCDFT PRIVATE
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    MKL::MKL
)

