cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
if(CMAKE_VERSION VERSION_LESS "4.0")
    cmake_policy(SET CMP0060 OLD)
endif()
cmake_policy(SET CMP0074 NEW)
project(steam C CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/cmake")

set(FORCE_STATIC_DEPS 0 CACHE BOOL "Force static linking of deps")
set(ENABLE_CUDA 0 CACHE BOOL "Enable CUDA")

if(FORCE_STATIC_DEPS)
    if(ENABLE_CUDA)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so CACHE INTERNAL "" FORCE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a CACHE INTERNAL "" FORCE)
        set(CMAKE_LINK_SEARCH_START_STATIC ON CACHE INTERNAL "" FORCE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()
    set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "" FORCE)
endif()

if(NOT XCODE AND NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Apply STEAM patches to the bundled mmseqs2 submodule before configuring it.
# Each patch lives in patches/*.patch and is applied at the lib/mmseqs root.
# Idempotent: if a patch is already applied (reverse-check succeeds) we skip
# it. If a patch cannot be applied and is also not already applied, the
# build aborts loudly — this happens when the mmseqs submodule pointer was
# bumped and the patch needs to be rebased against the new upstream code.
file(GLOB STEAM_PATCHES "${CMAKE_CURRENT_SOURCE_DIR}/patches/*.patch")
list(SORT STEAM_PATCHES)
find_package(Git QUIET REQUIRED)
foreach(PATCH ${STEAM_PATCHES})
    get_filename_component(PATCH_NAME ${PATCH} NAME)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} apply --reverse --check ${PATCH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs
        RESULT_VARIABLE _patch_reverse_ok
        OUTPUT_QUIET ERROR_QUIET
    )
    if(_patch_reverse_ok EQUAL 0)
        message(STATUS "STEAM patch already applied: ${PATCH_NAME}")
    else()
        execute_process(
            COMMAND ${GIT_EXECUTABLE} apply ${PATCH}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs
            RESULT_VARIABLE _patch_apply_ok
        )
        if(NOT _patch_apply_ok EQUAL 0)
            message(FATAL_ERROR
                "Failed to apply ${PATCH_NAME} to lib/mmseqs. "
                "If the mmseqs submodule pointer was bumped, the patch may "
                "need to be rebased against the new upstream code.")
        endif()
        message(STATUS "Applied STEAM patch: ${PATCH_NAME}")
    endif()
endforeach()

# Build MMseqs2 as a framework library only (no executable)
set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)

add_subdirectory(data)
add_subdirectory(src)
