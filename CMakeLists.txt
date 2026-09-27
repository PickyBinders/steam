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
# Idempotent: reverse the complete ordered series in a temporary source copy.
# Later patches may change earlier hunks, so independent reverse checks alone
# are insufficient. If the series is not applied and a patch cannot apply, the
# build aborts loudly — this happens when the mmseqs submodule pointer was
# bumped and the patch needs to be rebased against the new upstream code.
option(
    STEAM_SKIP_BUNDLED_PATCH_APPLICATION
    "Skip patch checks only when configuring an already-patched MMseqs worktree"
    OFF
)
file(GLOB STEAM_PATCHES "${CMAKE_CURRENT_SOURCE_DIR}/patches/*.patch")
list(SORT STEAM_PATCHES)
find_package(Git QUIET REQUIRED)
if(NOT STEAM_SKIP_BUNDLED_PATCH_APPLICATION)
    set(_patch_check_dir "${CMAKE_BINARY_DIR}/steam-patch-check")
    foreach(PATCH ${STEAM_PATCHES})
        file(READ "${PATCH}" _patch_text)
        string(REGEX MATCHALL "diff --git a/[^ \n]+" _patch_paths "${_patch_text}")
        foreach(_path ${_patch_paths})
            string(REPLACE "diff --git a/" "" _path "${_path}")
            get_filename_component(_parent "${_path}" DIRECTORY)
            file(MAKE_DIRECTORY "${_patch_check_dir}/${_parent}")
            configure_file("${CMAKE_CURRENT_SOURCE_DIR}/lib/mmseqs/${_path}"
                           "${_patch_check_dir}/${_path}" COPYONLY)
        endforeach()
    endforeach()
    set(_reverse_patches ${STEAM_PATCHES})
    list(REVERSE _reverse_patches)
    set(_series_applied TRUE)
    foreach(PATCH ${_reverse_patches})
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E env "GIT_CEILING_DIRECTORIES=${CMAKE_BINARY_DIR}"
                    ${GIT_EXECUTABLE} apply --reverse "${PATCH}"
            WORKING_DIRECTORY "${_patch_check_dir}"
            RESULT_VARIABLE _reverse_result OUTPUT_QUIET ERROR_QUIET
        )
        if(NOT _reverse_result EQUAL 0)
            set(_series_applied FALSE)
            break()
        endif()
    endforeach()
    if(_series_applied)
        message(STATUS "STEAM patch series already applied")
    else()
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
    endif()
endif()

# Build MMseqs2 as a framework library only (no executable)
set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE)
include(MMseqsSetupDerivedTarget)
add_subdirectory(lib/mmseqs EXCLUDE_FROM_ALL)

add_subdirectory(data)
add_subdirectory(src)
