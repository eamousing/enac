cmake_minimum_required(VERSION 3.2)

project(ENAC)
enable_language(Fortran)

set(sources 
    src/enac_commons.F90
    src/enac_functions.F90
    src/enac_main.F90
)

add_executable(enac ${sources})