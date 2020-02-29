### Initialize CMake (some internal stuff)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

project(NBI-TORQUE)
enable_language(Fortran)


### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)
include_directories(/usr/include)

message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("PROJLIBS = ${PROJLIBS}")

### Define flags for Fortran/C compiler
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(MY_DEBUG_FLAG_LIST
" -implicitnone "
" -g "
" -traceback "
#" -mkl "
" -O2 "
#  "-ggdb -C -p -pg -fbacktrace " 
#  "-ffpe-trap=invalid,zero,overflow,underflow " 
#  "-Werror " 
#  "-fmem-report " 
#  "-Waliasing -Wampersand  -Wline-truncation " 
#  "-Wnonstd-intrinsics  -Wsurprising -Wno-tabs -Wunderflow " 
#  "-Wall " 
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
find_package(SuiteSparse REQUIRED)
find_package(SuperLU REQUIRED)


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
	${SUPERLU_LIBRARIES}
	${SUITESPARSE_LIBRARIES}
	lapack
	#${lapack_lib}
	#${open_blas_lapack_lib}
	blas
	#${f77blas_lib} ${atlas_lib}
	#${open_blas_lib}
	debug efence
        )
