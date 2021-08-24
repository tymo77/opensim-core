set(EXAMPLE_2D_JUMPING_LINKAGE_FILES
        2D_linkage_leg.osim)

OpenSimAddExampleCXX(NAME example2DJumpingLinkage SUBDIR Moco
    EXECUTABLES example2DJumpingLinkage
    RESOURCES "${EXAMPLE_2D_JUMPING_LINKAGE_FILES}")

file(COPY ${EXAMPLE_2D_JUMPING_LINKAGE_FILES}
        DESTINATION
        "${CMAKE_BINARY_DIR}/Bindings/Java/Matlab/examples/Moco/example2DJumpingLinkage")
install(FILES ${EXAMPLE_2D_JUMPING_LINKAGE_FILES}
    DESTINATION "${OPENSIM_INSTALL_MATLABEXDIR}/Moco/example2DJumpingLinkage")
