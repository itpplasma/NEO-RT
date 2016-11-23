### Initialize CMake (some internal stuff)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(NBI-TORQUE)
enable_language(Fortran)

message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)
include_directories(/usr/include)

### Define flags for Fortran/C compiler
set(CMAKE_Fortran_FLAGS_RELEASE " -fimplicit-none")
set(MY_DEBUG_FLAG_LIST
" -fimplicit-none "
" -g "
#  "-ggdb -C -p -pg -fbacktrace " 
#  "-ffpe-trap=invalid,zero,overflow,underflow " 
#  "-Werror " 
#  "-fmem-report " 
#  "-Waliasing -Wampersand  -Wline-truncation " 
#  "-Wnonstd-intrinsics  -Wsurprising -Wno-tabs -Wunderflow " 
  "-Wall " 
#  "-Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation " 
  )
string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}") 
set(CFLAGS "-O2 -DDOUBLE_APPEND_FORTRAN")

### Define container for *.o and *.mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

### Build name extension
IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
  message("debug mode")
  set(EXE_NAME "nbi_torque.x.DEBUG")  
ELSE()
  message("release mode")
  set(EXE_NAME "nbi_torque.x")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)


### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(SuperLU_lib superlu_4.1 ${SUPERLU_DIR}lib NO_DEFAULT_PATH)
find_library(Umfpack_lib umfpack ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Amd_lib amd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Cholmod_lib cholmod ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Colamd_lib colamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Camd_lib camd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Metis_lib metis ${SUITESPARSE_DIR}lib ${METIS_PATH} NO_DEFAULT_PATH)
find_library(Ccolamd_lib ccolamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(lapack_lib lapack ${LAPACK_LIB} NO_DEFAULT_PATH)
find_library(f77blas_lib f77blas ${F77BLAS_LIB} NO_DEFAULT_PATH)
find_library(atlas_lib atlas ${ATLASBLAS_LIB} NO_DEFAULT_PATH)
find_library(open_blas_lapack_lib lapack ${OPEN_BLAS_LAPACK_LIB} NO_DEFAULT_PATH)
find_library(open_blas_lib blas ${OPEN_BLAS_LIB} NO_DEFAULT_PATH)

### Source files (load external file)
include (CMakeSources.txt)

### Define executable
add_executable(${EXE_NAME}
	${HYDRA_SRC_FILES}  
	${NEO_MAGFIE_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
)

### Add libraries linked with executable
target_link_libraries(${EXE_NAME} 
        ${SuperLU_lib}
        ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
	lapack
	#${lapack_lib}
	#${open_blas_lapack_lib}
	blas
	#${f77blas_lib} ${atlas_lib}
	#${open_blas_lib}
	debug efence
        )
