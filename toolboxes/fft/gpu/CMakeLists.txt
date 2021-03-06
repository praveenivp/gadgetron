if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUFFT__)
endif ()


cuda_add_library(gadgetron_toolbox_gpufft SHARED 
    cuNDFFT.h
    cuNDFFT.cpp
    cuNDFFT.cu
  )

set_target_properties(gadgetron_toolbox_gpufft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(gadgetron_toolbox_gpufft PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(gadgetron_toolbox_gpufft
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  gadgetron_toolbox_gpucore
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpufft
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  gpufft_export.h
  cuNDFFT.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)