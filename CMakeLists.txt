cmake_minimum_required(VERSION 2.6)
project(math_h)
option(tests "Compile tests.")
if(NOT DEFINED MATH_H_INC)
    get_directory_property(has_parent PARENT_DIRECTORY)
    if(has_parent)
	set(MATH_H_INC 
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
	    PARENT_SCOPE)
    endif()

    include_directories(include)

    file(GLOB EXAMPLES Example/*.cpp)
    foreach(example ${EXAMPLES})
	get_filename_component(ex_app ${example} NAME_WE)
	add_executable(${ex_app} ${example})
    endforeach()

    if(tests)
	file(GLOB TESTS tests/*.cpp)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
	find_package(GTest REQUIRED)
	include_directories(${GTEST_INCLUDE_DIRS})
	add_executable(math_h_test.exe ${TESTS})
	target_link_libraries(math_h_test.exe ${GTEST_BOTH_LIBRARIES})
    endif()
endif()
