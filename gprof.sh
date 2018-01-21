#!/usr/bin/env bash

set -ex

# grep CONFIG_HZ= /boot/config*

PROJECT_ROOT="$(dirname $0)"

PROFILEFREQUENCY=1000 LD_PRELOAD=/usr/lib/libprofiler.so.0 CPUPROFILE=/tmp/sweepline.profile "$PROJECT_ROOT/sweepline"
google-pprof --pdf "$PROJECT_ROOT/sweepline" /tmp/sweepline.profile > /tmp/sweepline.pdf
gnome-open /tmp/sweepline.pdf
