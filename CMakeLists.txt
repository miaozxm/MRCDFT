cmake_minimum_required(VERSION 3.15)

project(MRCDFT LANGUAGES Fortran)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(EXE_DIR "${CMAKE_SOURCE_DIR}/bin")
set(MOD_DIR "${CMAKE_BINARY_DIR}/mod")
file(GLOB SRC_FILES CONFIGURE_DEPENDS "${SRC_DIR}/*.f90")


# Create the executable target
add_executable(MRCDFT ${SRC_FILES})
set_target_properties(MRCDFT PROPERTIES
    Fortran_MODULE_DIRECTORY "${MOD_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${EXE_DIR}"
)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(CompilerChecks)

include(CompilerOptions)
target_compile_options(MRCDFT PRIVATE ${COMPILER_OPTIONS})
target_link_options(MRCDFT PRIVATE ${LINKER_OPTIONS})

include(Dependencies)
target_link_libraries(MRCDFT PRIVATE
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    MKL::MKL
)

include(Installation)
include(Packaging)