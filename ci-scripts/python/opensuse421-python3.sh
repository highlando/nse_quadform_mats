#!/bin/bash

# install packages via package manager
set -ex
zypper -n install python3-devel python3-numpy python3-scipy


# execute python scripts
cd python


# cylinderwake steadystate
python cylinderwake_steadystate.py
python cylinderwake_steadystate.py --N=1 --Re=50 --Picardsteps=5
python cylinderwake_steadystate.py --N=2 --Re=50 --Picardsteps=5
#python cylinderwake_steadystate.py --N=3 --Re=50 --Picardsteps=5


# cylinderwake steadystate bccontrol
python cylinderwake_steadystate_bccontrol.py
python cylinderwake_steadystate_bccontrol.py --N=1 --Re=50 --Picardsteps=5 --palpha=1e-4
python cylinderwake_steadystate_bccontrol.py --N=2 --Re=50 --Picardsteps=5 --palpha=1e-4
#python cylinderwake_steadystate_bccontrol.py --N=3 --Re=50 --Picardsteps=5 --palpha=1e-4


# cylinderwake tdp pout vout
python cylinderwake_tdp_pout_vout.py --N=1 --Re=50 --Picardsteps=5 --t0=0.0 --tE=1.5 --Nts=4096
python cylinderwake_tdp_pout_vout.py --N=2 --Re=50 --Picardsteps=5 --t0=0.0 --tE=1.5 --Nts=4096
#python cylinderwake_tdp_pout_vout.py --N=3 --Re=50 --Picardsteps=5 --t0=0.0 --tE=1.5 --Nts=4096


# cylinderwake tdp pout vout bccontrol
python cylinderwake_tdp_pout_vout_bccontrol.py
python cylinderwake_tdp_pout_vout_bccontrol.py --N=1 --Re=50 --Picardsteps=5 --palpha=1e-4 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096
python cylinderwake_tdp_pout_vout_bccontrol.py --N=2 --Re=50 --Picardsteps=5 --palpha=1e-4 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096
#python cylinderwake_tdp_pout_vout_bccontrol.py --N=3 --Re=50 --Picardsteps=5 --palpha=1e-4 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096


# drivencavity steadystate
python drivencavity_steadystate.py
python drivencavity_steadystate.py --N=10 --Re=500 --Picardsteps=5
python drivencavity_steadystate.py --N=20 --Re=500 --Picardsteps=5
python drivencavity_steadystate.py --N=30 --Re=500 --Picardsteps=5


# drivencavity tdp pout vout distrcontrol
python drivencavity_tdp_pout_vout_distcontrol.py 
python drivencavity_tdp_pout_vout_distcontrol.py --N=10 --Re=500 --Picardsteps=5 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096
python drivencavity_tdp_pout_vout_distcontrol.py --N=20 --Re=500 --Picardsteps=5 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096
python drivencavity_tdp_pout_vout_distcontrol.py --N=30 --Re=500 --Picardsteps=5 --omega=3.5 --t0=0.0 --tE=1.5 --Nts=4096


# cylinderwake stokes
python stokes_ss.py
python stokes_ss.py --N=1 --RE=50 --Picardsteps=5 
python stokes_ss.py --N=2 --RE=50 --Picardsteps=5 
#python stokes_ss.py --N=3 --RE=50 --Picardsteps=5 

# test units
python test_units.py

