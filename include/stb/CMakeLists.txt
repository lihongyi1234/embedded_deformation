cmake_minimum_required(VERSION 3.1)

# Create a library for the viewer code
add_library(
    stb 
    stb_impl.cpp
    )

target_include_directories(stb PRIVATE "../include")    

set_target_properties(stb PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
