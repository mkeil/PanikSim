cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

find_package(GLog REQUIRED)
include_directories(${GLOG_INCLUDE_DIR})

find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})


include_directories("src")

set (CUDA_NVCC_FLAGS "-arch=sm_20")
	
cuda_add_executable(PanicSim src/PanicSimMain.cu src/hostFunc.cu src/kernels.cu src/Xlib_mod.cu) 
target_link_libraries(PanicSim ${X11_LIBRARIES})
target_link_libraries(PanicSim ${GLOG_LIBRARY})

cuda_add_executable(PanicSim_constant_stepsize src/PanicSimMain.cu src/hostFunc.cu src/kernels_constant_stepsize.cu src/Xlib_mod.cu) 
target_link_libraries(PanicSim_constant_stepsize ${X11_LIBRARIES})
target_link_libraries(PanicSim_constant_stepsize ${GLOG_LIBRARY})

# common flags
SET(CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-long-long -Wno-variadic-macros")
SET(CMAKE_CXX_FLAGS_RELEASE "-DNEXPENSIVE -O2 -g0")
SET(CMAKE_CXX_FLAGS_DEBUG "-DWARNINGS -O0 -g3 -fno-inline")