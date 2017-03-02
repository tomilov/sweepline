project("sweepline")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

#[[
if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "only clang supported currently")
endif()

include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
if(NOT COMPILER_SUPPORTS_LIBCXX)
    message(FATAL_ERROR "libc++ is not installed or libc++ not supported by clang++")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi -lunwind")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # suppress -rdynamic flag

add_compile_options(-fconstexpr-steps=50000000 -ftemplate-depth=200)
add_compile_options(-ftemplate-backtrace-limit=0 -fconstexpr-backtrace-limit=0 -fdiagnostics-show-template-tree)
add_compile_options(-W -Weverything -Wconversion -Wmismatched-tags -Wmissing-include-dirs -Wno-c++98-compat -Wno-c++98-compat-pedantic)
add_compile_options(-W -Wconversion -Wmissing-include-dirs -Wno-c++98-compat -Wno-c++98-compat-pedantic)
add_compile_options(-Wno-gnu-anonymous-struct -Wno-nested-anon-types)
add_compile_options(-Wno-c++1z-extensions)
add_compile_options(-Wuninitialized -Wfatal-errors -Wno-padded)
#add_compile_options(-ftime-report)

#]]

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -fno-inline -DDEBUG=1 -D_GLIBCXX_DEBUG=1 -march=x86-64 -mtune=generic")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -gline-tables-only -fno-omit-frame-pointer -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG=1 -march=native")
else()
    message(STATUS "Wrong build type selected, default to Debug.")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -gline-tables-only -fno-omit-frame-pointer")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")

include_directories(".")

add_executable("sweepline" "main.cpp" "sweepline.hpp" "rb_tree.hpp")

set_target_properties("sweepline" PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
)
