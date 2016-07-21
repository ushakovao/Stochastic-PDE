cmake_minimum_required(VERSION 2.8)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
        find_package(Feel++
                PATHS $ENV{FEELPP_DIR}/share/feel/cmake/modules
                /usr/share/feel/cmake/modules
                /usr/local/share/feel/cmake§modules
                opt/share/feel/cmake/modules
        )
        if(NOT FEELPP_FOUND)
                message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the illation directory.")
        endif()
endif()
feelpp_add_application(
 BS1d
  SRCS BlackScholes1d.cpp
  GEO square1d.geo
  DEFS A_DEF=2
  CFG BlackScholes1d.cfg )

