cmake_minimum_required(VERSION 2.8)

include(PICMake)

set(CMAKE_CXX_STANDARD 11)

INCLUDE_DIRECTORIES(src src/exiv2 src/xmpsdk/include)

pi_add_target(svar_exif SHARED src REQUIRED System DEPENDENCY -lz)

pi_install(TARGETS svar_exif)
