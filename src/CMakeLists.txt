add_library( myslam SHARED
    frame.cpp
    camera.cpp
    config.cpp
    g2o_types.cpp
    simulator.cpp
    queueSimulator.cpp
)

target_link_libraries( myslam
    ${THIRD_PARTY_LIBS}
)
