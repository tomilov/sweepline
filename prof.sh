#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/.

BUILD_DIR=$PROJECT_ROOT/.
cat /tmp/s | PROFILEFREQUENCY=10000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=$BUILD_DIR/sweepline.profile $BUILD_DIR/sweepline
google-pprof --pdf $BUILD_DIR/sweepline $BUILD_DIR/sweepline.profile > $BUILD_DIR/prof.pdf
gnome-open $BUILD_DIR/prof.pdf
