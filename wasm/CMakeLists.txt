cmake_minimum_required(VERSION 3.10)
project(chapro_wasm C)

set(CHAPRO_SOURCES
    ../cha_core.c ../cha_scale.c ../db.c ../fft.c ../rfft.c
    ../agc_prepare.c ../agc_process.c
    ../cfirfb_prepare.c ../cfirfb_process.c
    ../firfb_prepare.c ../firfb_process.c
    ../iirfb_design.c ../iirfb_prepare.c ../iirfb_process.c
    ../ciirfb_design.c ../ciirfb_prepare.c ../ciirfb_process.c
    ../afc_prepare.c ../afc_filters.c ../afc_process.c
    ../nfc_prepare.c ../nfc_process.c
    ../sha_prepare.c ../sha_process.c
    ../icmp_prepare.c ../icmp_process.c
)

add_executable(chapro_wasm ${CHAPRO_SOURCES} main.c)
set_target_properties(chapro_wasm PROPERTIES OUTPUT_NAME "chapro")

target_include_directories(chapro_wasm PRIVATE ../)

# Emscripten specific link options
target_link_options(chapro_wasm PRIVATE
    "-s" "EXPORTED_FUNCTIONS=['_cha_version']"
    "-s" "EXPORTED_RUNTIME_METHODS=['ccall','cwrap']"
    "-s" "MODULARIZE=1"
    "-s" "EXPORT_NAME='createChapro'"
)
