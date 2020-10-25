# CmakeLists.txt in CrossMIPs dir

# Adds a library called crossmips (crossmips.a under Linux, crossmips.lib under Windows) from the all .h and .cpp files
file(GLOB crossmips_headers *.h)
file(GLOB crossmips_sources *.cpp)

# check CUDA option
option (WITH_CUDA "Support to CUDA" OFF) 
if(WITH_CUDA)

	message(STATUS "Checking CUDA")
	find_package(CUDA)
	
	if(CUDA_FOUND)
	
		# enable CUDA-dependent code
		add_definitions(-DUSECUDA)

		CUDA_ADD_LIBRARY(crossmips STATIC ${crossmips_headers} ${crossmips_sources} compute_funcs.cu OPTIONS -arch sm_30)
		
	else(CUDA_FOUND)
	
		message(FATAL_ERROR "Checking CUDA ... error!!!")
		
	endif(CUDA_FOUND)
	
else(WITH_CUDA)

	add_library(crossmips STATIC ${crossmips_headers} ${crossmips_sources})

endif(WITH_CUDA)

set_target_properties(crossmips PROPERTIES POSITION_INDEPENDENT_CODE ON)

TARGET_LINK_LIBRARIES(crossmips ${CUDA_LIBRARIES})
