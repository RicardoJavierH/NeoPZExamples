set(tutorial_sources "")

add_subdirectory(Poisson2D)
add_subdirectory(Projections)
add_subdirectory(HCurl3D)
add_subdirectory(ExperimentosLCC)


#enable folder structure for the tutorials (needed in some IDEs)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${tutorial_sources})