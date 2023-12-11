cmake_minimum_required(VERSION 3.21)
project(analysis_k40)

set(CMAKE_CXX_STANDARD 17)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree)

# Include ROOT's headers
include(${ROOT_USE_FILE})

# Specify the ROOT headers directory
#include_directories(/lustre/neutrino/huangweilun/ROOT/root/include)

add_executable(analysis analysis_trigger_root_tts.cpp analysis_trigger_root_tts.h)
target_link_libraries(analysis ${ROOT_LIBRARIES})
