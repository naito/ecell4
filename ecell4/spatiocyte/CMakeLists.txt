if (NO_SHARED)
    return()
endif()

file(GLOB CPP_FILES *.cpp)
file(GLOB HPP_FILES *.hpp)

add_library(ecell4-spatiocyte SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-spatiocyte ecell4-core)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-spatiocyte DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/spatiocyte")
