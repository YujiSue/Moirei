cmake_minimum_required(VERSION 3.14)
project(moirei)
include(ExternalProject)
if(MSVC)
	add_compile_definitions(_USE_DLL)
	set(CMAKE_CXX_FLAGS "/utf-8")
	link_libraries(sobjd.lib sbioinfod.lib sappd.lib)
elseif(UNIX)
	if (DEFINED ENV{SLIB_INCLUDE})
		message(STATUS "SLIB_INCLUDE=$ENV{SLIB_INCLUDE}")
	else()
		#message(STATUS "SLIB_INCLUDE not defined.")
		set(ENV{SLIB_INCLUDE} ${CMAKE_INSTALL_PREFIX}/include/slib)
		message(STATUS "SLIB_INCLUDE=$ENV{SLIB_INCLUDE}")
	endif()
	if (DEFINED ENV{SLIB_DIR})
		message(STATUS "SLIB_DIR=$ENV{SLIB_DIR}")
	else()
		#message(STATUS "SLIB_DIR not defined.")
		set(ENV{SLIB_DIR} ${CMAKE_INSTALL_PREFIX}/lib)
		message(STATUS "SLIB_DIR=$ENV{SLIB_DIR}")		
	endif()
	find_library(LSOBJ_INSTALLED NAMES sobj)
	find_library(LSBI_INSTALLED NAMES sbioinfo)
	find_library(LSAPP_INSTALLED NAMES sapp)
	if (LSOBJ_INSTALLED AND LSBI_INSTALLED AND LSAPP_INSTALLED)
		message(STATUS "slib found.")
	else()
		message(STATUS "Error: Three libraries <libsobj, libsbioinfo, libsapp> are required.")
		ExternalProject_Add(
			slib
			PREFIX ${CMAKE_BINARY_DIR}/ext
			GIT_REPOSITORY https://github.com/YujiSue/slib.git
			GIT_TAG master
    			CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/ext -DBUILD_BIOINFO_LIB=ON -DBUILD_SAPP_LIB=ON 
    			UPDATE_DISCONNECTED 1
		)
		ExternalProject_Get_Property(slib install_dir)
		include_directories(${install_dir}/include)
		link_directories(${install_dir}/lib)	
	endif()
	message(STATUS "Install to ${CMAKE_INSTALL_PREFIX}")
	set(CMAKE_CXX_FLAGS "-std=c++11 -O2")
	link_libraries(sobj sbioinfo sapp dl)
endif()
include_directories(${PROJECT_SOURCE_DIR}/src $ENV{SLIB_INCLUDE})
link_directories($ENV{SLIB_DIR})
execute_process(
	COMMAND python MakeProf.py profile.h.prv appinfo.json ./src/profile.h
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_subdirectory(src)
