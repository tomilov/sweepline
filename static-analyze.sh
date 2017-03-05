#!/usr/bin/env bash

set -vex

#git clean -fdX
mkdir -p build/
pushd build/
cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_COLOR_MAKEFILE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_VERBOSE_MAKEFILE=ON ..
cmake --build .
popd

#make --jobs $((`nproc` + 1))

clang-tidy -p build/ main.cpp -header-filter=.* -checks=*,-llvm-header-guard,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-readability-else-after-return,-clang-analyzer-alpha.core.PointerArithm,-clang-analyzer-alpha.clone.CloneChecker,-readability-implicit-bool-cast,-llvm-namespace-comment,-llvm-include-order,-google-readability-namespace-comments,-google-runtime-references,-cppcoreguidelines-pro-type-union-access
