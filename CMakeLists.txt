cmake_minimum_required(VERSION 3.7)
project(particles)

set(CMAKE_C_STANDARD 11)


set(SOURCE_FILES
        main.c
        common_use.c
        common_use.h
        dust_wave.c
        dust_wave.h
        gas_wave.c
        gas_wave.h
        explicit.c
        explicit.h
        xy_system.c
        xy_system.h
        cells.c
        cells.h)
add_executable(particles ${SOURCE_FILES})
target_link_libraries(particles m)