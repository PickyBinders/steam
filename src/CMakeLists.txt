include_directories(commons)
include_directories(.)
add_subdirectory(alignment)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(workflow)

add_library(steam-framework
        ${alignment_source_files}
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        SteamBase.cpp
)
mmseqs_setup_derived_target(steam-framework)
target_include_directories(steam-framework PUBLIC ${CMAKE_BINARY_DIR}/generated)
add_dependencies(steam-framework local-generated)

add_subdirectory(version)
add_executable(steam steam.cpp)
mmseqs_setup_derived_target(steam steam-framework)
target_link_libraries(mmseqs-framework version)
target_link_libraries(steam version)
install(TARGETS steam DESTINATION bin)
