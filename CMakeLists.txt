cmake_minimum_required (VERSION 3.1)

project (ECAL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})
include(${ROOT_USE_FILE})

add_executable(calibrate "calibrate.cpp")
target_link_libraries(calibrate ${ROOT_LIBRARIES})

set(SHELL_SCRIPTS
    dat_to_root.sh
    mk_dir.sh
    txt_list.sh
    ${pyscripts}
)

foreach(script ${SHELL_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${script}
        ${PROJECT_BINARY_DIR}/${script}
        COPYONLY
    )
endforeach()


