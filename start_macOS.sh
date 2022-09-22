#!/bin/bash
#------------------------------------------------------------------
# set relative path from run.sh
#------------------------------------------------------------------
path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd ${path}
#------------------------------------------------------------------
# execute .jl code(s)
#------------------------------------------------------------------
#julia -i -O3 --threads=1 --project=.
julia -i -O3 --threads=1 --check-bounds=no --project=.
#julia -i -O3 -t auto --check-bounds=no --project=.