#!/bin/bash

# install packages via package manager
set -e 
apt-get update  --yes
apt-get install --yes python-dev python-scipy

# execute python scripts
cd python


# cylinderwake steadystate
python cylinderwake_steadystate.py

# cylinderwake steadystate bccontrol
python cylinderwake_steadystate_bccontrol.py

# cylinderwake tdp pout vout
python cylinderwake_tdp_pout_vout.py

# cylinderwake tdp pout vout bccontrol
python cylinderwake_tdp_pout_vout_bccontrol.py

# cylinderwake stokes
python stokes_ss.py

# drivencavity steadystate
python drivencavity_steadystate.py

# drivencavity tdp pout vout distrcontrol
python drivencavity_tdp_pout_vout_distcontrol.py

# test units
python test_units.py
