cmake_minimum_required(VERSION 3.10)

project(sparse)

option(Enable_BLAS "Enable BLAS" OFF)

file(GLOB source_files
            "*.c"
            "src/*.c"
            "include/*.h"
    )


include_directories(include)
add_executable(sparse ${source_files})
target_link_libraries(sparse PUBLIC m)

if(Enable_BLAS)
    find_package(BLAS)
    if(BLAS_FOUND)
        target_link_libraries(sparse PRIVATE ${BLAS_LIBRARIES})
        SET(HAVE_CBLAS 1)
    endif()
endif()
include_directories(include)

configure_file(config.h.in ${CMAKE_SOURCE_DIR}/include/config.h @ONLY)






