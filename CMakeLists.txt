cmake_minimum_required(VERSION 3.10)
project(chapro C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Source files for the CHAPRO library
set(CHAPRO_SOURCES
    cha_core.c cha_scale.c db.c fft.c rfft.c
    agc_prepare.c agc_process.c
    cfirfb_prepare.c cfirfb_process.c
    firfb_prepare.c firfb_process.c
    iirfb_design.c iirfb_prepare.c iirfb_process.c
    ciirfb_design.c ciirfb_prepare.c ciirfb_process.c
    afc_prepare.c afc_filters.c afc_process.c
    nfc_prepare.c nfc_process.c
    sha_prepare.c sha_process.c
    icmp_prepare.c icmp_process.c
)

# Build both static and shared libraries
add_library(chapro_static STATIC ${CHAPRO_SOURCES})
add_library(chapro SHARED ${CHAPRO_SOURCES})

set_target_properties(chapro_static PROPERTIES OUTPUT_NAME chapro)

# Find dependencies: math, zlib, sigpro, arsc
find_library(M_LIB m)
find_library(Z_LIB z)
find_library(SIGPRO_LIB sigpro HINTS /usr/local/lib)
find_library(ARSC_LIB arsc HINTS /usr/local/lib)

# Try to find headers
find_path(SIGPRO_INCLUDE_DIR sigpro.h HINTS /usr/local/include)
find_path(ARSC_INCLUDE_DIR arsclib.h HINTS /usr/local/include)

# Include directories
target_include_directories(chapro_static PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(chapro PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

if (SIGPRO_INCLUDE_DIR)
    target_include_directories(chapro_static PRIVATE ${SIGPRO_INCLUDE_DIR})
    target_include_directories(chapro PRIVATE ${SIGPRO_INCLUDE_DIR})
endif()

if (ARSC_INCLUDE_DIR)
    target_include_directories(chapro_static PRIVATE ${ARSC_INCLUDE_DIR})
    target_include_directories(chapro PRIVATE ${ARSC_INCLUDE_DIR})
endif()

# Dependencies
set(CHAPRO_DEPS ${M_LIB})
if (Z_LIB)
    list(APPEND CHAPRO_DEPS ${Z_LIB})
endif()
if (SIGPRO_LIB)
    list(APPEND CHAPRO_DEPS ${SIGPRO_LIB})
endif()
if (ARSC_LIB)
    list(APPEND CHAPRO_DEPS ${ARSC_LIB})
endif()

target_link_libraries(chapro PRIVATE ${CHAPRO_DEPS})

# Testing
include(CTest)
if(BUILD_TESTING)
    set(TEST_PGMS
        tst_cffa tst_cffio tst_cffsc tst_ffa tst_ffio tst_ffsc
        tst_cifa tst_cifio tst_cifsc tst_ifa tst_ifio tst_ifsc
        tst_gha tst_nfc tst_sha gha_demo opt_afc
    )

    foreach(test_name ${TEST_PGMS})
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.c")
            add_executable(${test_name} ${test_name}.c)
            target_link_libraries(${test_name} PRIVATE chapro_static ${CHAPRO_DEPS})
            if(APPLE)
                target_link_libraries(${test_name} PRIVATE "-framework CoreServices" "-framework CoreAudio" "-framework AudioUnit" "-framework AudioToolbox")
            elseif(UNIX AND NOT APPLE)
                target_link_libraries(${test_name} PRIVATE asound)
            endif()
            add_test(NAME ${test_name} COMMAND ${test_name})
        endif()
    endforeach()
endif()
