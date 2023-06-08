if(NOT dune-surface_int_FOUND)
# Whether this module is installed or not
set(dune-surface_int_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/zavala68/mydune/dune-surface_int)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-surface_int_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-surface_int_INCLUDE_DIRS "/home/zavala68/mydune/dune-surface_int")
set(dune-surface_int_CXX_FLAGS "-std=c++17 ")
set(dune-surface_int_CXX_FLAGS_DEBUG "-g")
set(dune-surface_int_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-surface_int_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-surface_int_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-surface_int_DEPENDS "dune-common;dune-uggrid;dune-geometry;dune-localfunctions;dune-grid;dune-istl;dune-typetree;dune-functions;dune-foamgrid;dune-vtk;dune-curvedgeometry;dune-curvedgrid;dune-gmsh4")
set(dune-surface_int_SUGGESTS "")
set(dune-surface_int_MODULE_PATH "/home/zavala68/mydune/dune-surface_int/cmake/modules")
set(dune-surface_int_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-surface_int_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-surface_int-targets.cmake")
endif()
endif()
