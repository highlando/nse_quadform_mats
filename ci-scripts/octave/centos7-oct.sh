#!/bin/bash

# install packages via package manager
set -ex 
yum makecache fast
yum install -y octave


# execute matlab scripts
cd matlab


# cylinderwake steadystate
octave-cli  --eval "add_to_path; cylinderwake_steadystate"
octave-cli  --eval "add_to_path; cylinderwake_steadystate(1, 50, 5)"
octave-cli  --eval "add_to_path; cylinderwake_steadystate(2, 50, 5)"
octave-cli  --eval "add_to_path; cylinderwake_steadystate(3, 50, 5)"


# cylinderwake steadystate bccontrol
octave-cli  --eval "add_to_path; cylinderwake_steadystate_bccontrol"
octave-cli  --eval "add_to_path; cylinderwake_steadystate_bccontrol(1, 50, 5, 1e-4)"
octave-cli  --eval "add_to_path; cylinderwake_steadystate_bccontrol(2, 50, 5, 1e-4)"
octave-cli  --eval "add_to_path; cylinderwake_steadystate_bccontrol(3, 50, 5, 1e-4)"


# cylinderwake tdp pout vout
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout(1, 50, 5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout(2, 50, 5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout(3, 50, 5, 0.0, 1.5, 4096)"


# cylinderwake tdp pout vout bccontrol
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout_bccontrol"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout_bccontrol(1, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout_bccontrol(2, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; cylinderwake_tdp_pout_vout_bccontrol(3, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096)"


# drivencavity steadystate
octave-cli  --eval "add_to_path; drivencavity_steadystate"
octave-cli  --eval "add_to_path; drivencavity_steadystate(10, 500, 5)"
octave-cli  --eval "add_to_path; drivencavity_steadystate(20, 500, 5)"
octave-cli  --eval "add_to_path; drivencavity_steadystate(30, 500, 5)"


# drivencavity tdp pout vout distrcontrol
octave-cli  --eval "add_to_path; drivencavity_tdp_pout_vout_distcontrol"
octave-cli  --eval "add_to_path; drivencavity_tdp_pout_vout_distcontrol(10, 500, 5, 3.5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; drivencavity_tdp_pout_vout_distcontrol(20, 500, 5, 3.5, 0.0, 1.5, 4096)"
octave-cli  --eval "add_to_path; drivencavity_tdp_pout_vout_distcontrol(30, 500, 5, 3.5, 0.0, 1.5, 4096)"


# cylinderwake stokes
octave-cli  --eval "add_to_path; stokes_ss"
octave-cli  --eval "add_to_path; stokes_ss(1, 50, 5)"
octave-cli  --eval "add_to_path; stokes_ss(2, 50, 5)"
octave-cli  --eval "add_to_path; stokes_ss(3, 50, 5)"


