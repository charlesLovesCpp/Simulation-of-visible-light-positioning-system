add_executable( main main.cpp )
target_link_libraries( main myslam ${Sophus_LIBRARIES} )

#target_link_libraries( main myslam ${Sophus_LIBRARIES} ${CSPARSE_LIBRARY} g2o_core g2o_stuff g2o_types_sba g2o_csparse_extension)
#include_directories( ${PROJECT_SOURCE_DIR}/include )