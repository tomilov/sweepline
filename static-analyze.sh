#!/usr/bin/env bash

set -vex

PROJECT_ROOT="$( dirname "$( realpath "$0" )" )"

pushd "$( mktemp -d )"

cmake -D CMAKE_CXX_COMPILER="$( which clang++-5.0 )" -D CMAKE_COLOR_MAKEFILE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_VERBOSE_MAKEFILE=ON "$PROJECT_ROOT"

cmake --build . -- --jobs $(( $( nproc ) + 1 ))

clang-tidy-5.0 -p . "$PROJECT_ROOT/main.cpp" -header-filter=.* -checks=*,-llvm-header-guard,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-readability-else-after-return,-clang-analyzer-alpha.core.PointerArithm,-clang-analyzer-alpha.clone.CloneChecker,-readability-implicit-bool-cast,-llvm-namespace-comment,-llvm-include-order,-google-readability-namespace-comments,-google-runtime-references,-cppcoreguidelines-pro-type-union-access

popd
