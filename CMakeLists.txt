cmake_minimum_required(VERSION 2.8)

project(SCAFFOLD_GRAPH_SET)

include_directories(${SCAFFOLD_GRAPH_SOURCE_DIR})

add_executable(scaffold_graph
               main.cpp
               dotwriter.cpp
               dotwriter.h
               graph.cpp
               graph.h
               read.cpp
               read.h
               sam_reader.cpp
               sam_reader.h
               )
			   
target_link_libraries(scaffold_graph libyaml-cpp.a)			   

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -pedantic -O2 -std=c++11")
     # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic -std=c++1z -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2 -fno-sanitize-recover -fstack-protector")
endif()
