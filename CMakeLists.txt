cmake_minimum_required( VERSION 2.8 )
project( myslam )

set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )
# set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
# set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib )

############### dependencies ######################
# Eigen
include_directories( " /usr/include/eigen3" )

# OpenCV
find_package( OpenCV 3.1 REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

# sophus
find_package( Sophus REQUIRED )
include_directories( ${Sophus_INCLUDE_DIRS} )

# g2o
find_package( G2O REQUIRED )
include_directories( ${G2O_INCLUDE_DIRS} )

find_package( CSparse REQUIRED )
include_directories(${CSPARSE_INCLUDE_DIR})

find_package( Cholmod REQUIRED )
include_directories( ${CHOLMOD_INCLUDE_DIR} )

# Ceres 
# find_package( Ceres REQUIRED )
# include_directories( ${CERES_INCLUDE_DIRS} )

set( THIRD_PARTY_LIBS
      ${OpenCV_LIBS}
      ${CHOLMOD_LIBRARIES}
      ${Sophus_LIBARIES}
      ${CSPARSE_INCLUDE_DIR}
      ${CSPARSE_LIBRARY}
      g2o_core g2o_stuff g2o_types_sba g2o_csparse_extension)
############### dependencies ######################
include_directories( ${PROJECT_SOURCE_DIR}/include )
add_subdirectory( src )
add_subdirectory( test )

#add_executable( main main.cpp )
#target_link_libraries( main g2o_core g2o_stuff g2o_types_slam3d ${CHOLMOD_LIBRARIES} ${Sophus_LIBRARIES} ${OpenCV_LIBS})

