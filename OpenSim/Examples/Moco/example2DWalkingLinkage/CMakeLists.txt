set(EXAMPLE_2D_WALKING_LINKAGE_FILES
        2D_gait_linkage.osim
        referenceCoordinates.sto
        referenceGRF.sto
        referenceGRF.xml)

OpenSimAddExampleCXX(NAME example2DWalkingLinkage SUBDIR Moco
    EXECUTABLES example2DWalkingLinkage
    RESOURCES "${EXAMPLE_2D_WALKING_LINKAGE_FILES}")

file(COPY ${EXAMPLE_2D_WALKING_LINKAGE_FILES}
        DESTINATION
        "${CMAKE_BINARY_DIR}/Bindings/Java/Matlab/examples/Moco/example2DWalkingLinkage")
install(FILES ${EXAMPLE_2D_WALKING_LINKAGE_FILES}
    DESTINATION "${OPENSIM_INSTALL_MATLABEXDIR}/Moco/example2DWalkingLinkage")
