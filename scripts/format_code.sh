#!/usr/bin/env bash

module reset
module load compiler Clang
clang-format -i --style=Microsoft $(ls src/*/*.cpp src/*/*.hpp | egrep -v "src/host/xcl2\..*|src/test/catch_amalgamated\..*|src/host/json\.hpp")