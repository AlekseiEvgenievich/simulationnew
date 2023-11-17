#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CrystalEye)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
#SET(CADmesh_lib /home/DAMPE/libowu/local/lib)
		
SET(CADmesh_lib /Users/aleksejsmirnov/Software/CADmesh/install/lib)
SET(CADmesh_inc /Users/aleksejsmirnov/Software/CADmesh/install/include)
#	SET(CADmesh_inc /home/DAMPE/libowu/local/include)
include_directories(${CADmesh_inc})
link_directories(${CADmesh_lib})

#----------------------------------------------------------------------------


# CADMesh
find_package(cadmesh REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(
		${PROJECT_SOURCE_DIR}/include
		${Geant4_INCLUDE_DIR}
		${ROOT_INCLUDE_DIR}
		${CADMESH_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB macros RELATIVE ${PROJECT_SOURCE_DIR} macros/*.mac)
file(GLOB condor RELATIVE ${PROJECT_SOURCE_DIR} condor/*.sh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(CrystalEye CrystalEye.cc ${sources} ${headers})
target_link_libraries(CrystalEye ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${cadmesh_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B3. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

foreach(_script ${condor})
		configure_file(
				${_script}
				${PROJECT_BINARY_DIR}/${_script}
				COPYONLY
		)
endforeach()

foreach(_script submit_theta_phi.sh submit_bkg.sh submit_theta.sh)
		configure_file(
				${_script}
				${PROJECT_BINARY_DIR}/${_script}
				COPYONLY
		)
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
##add_custom_target(CrystalEye DEPENDS CrystalEye)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CrystalEye DESTINATION bin)
