project("sweepline")

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # suppress -rdynamic flag

    #[[
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
    if(NOT COMPILER_SUPPORTS_LIBCXX)
        message(FATAL_ERROR "libc++ is not installed or libc++ not supported by clang++")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    ]]

    #[[
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -gline-tables-only -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=undefined")
    ]]

    add_compile_options(-W -Weverything)
    add_compile_options(-fconstexpr-steps=50000000 -fconstexpr-backtrace-limit=0 -fdiagnostics-show-template-tree)
    add_compile_options(-Wmismatched-tags)
    add_compile_options(-Wno-c++98-compat -Wno-c++98-compat-pedantic)
    add_compile_options(-Wno-gnu-anonymous-struct -Wno-nested-anon-types)
    add_compile_options(-Wno-c++1z-extensions)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-W -Wall)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} currently is not supported")
endif()

add_compile_options(-ftemplate-depth=200)
add_compile_options(-ftemplate-backtrace-limit=0)
add_compile_options(-Wconversion -Wmissing-include-dirs)
add_compile_options(-Wuninitialized -Wfatal-errors -Wno-padded)

#add_compile_options(-ftime-report)

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -fno-inline -DDEBUG=1 -D_GLIBCXX_DEBUG=1 -march=x86-64 -mtune=generic")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-fno-omit-frame-pointer -DNDEBUG=1 -march=native")
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -gline-tables-only")
    elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Og -ggdb")
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG=1 -march=native")
else()
    message(STATUS "Wrong build type selected, defaulted to Debug.")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

include_directories(".")

add_executable("sweepline" "main.cpp" "sweepline.hpp" "rb_tree.hpp")

set_target_properties("sweepline" PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
    )
