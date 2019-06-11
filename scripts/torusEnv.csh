# Set up paths for Torus under C-like shells
# Call with: source scripts/torusEnv.csh
set this_dir=`pwd`
setenv PATH ${this_dir}/utilties/python:${this_dir}/bin:${this_dir}/scripts:${PATH}
