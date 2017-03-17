#!/bin/bash

# install packages via package manager
set -e 
apt-get update  --yes
apt-get install --yes python3-dev python3-scipy

# execute python scripts
cd python


# cylinderwake steadystate
python3 cylinderwake_steadystate.py

# cylinderwake steadystate bccontrol
python3 cylinderwake_steadystate_bccontrol.py

# cylinderwake tdp pout vout
python3 cylinderwake_tdp_pout_vout.py

# cylinderwake tdp pout vout bccontrol
python3 cylinderwake_tdp_pout_vout_bccontrol.py

# cylinderwake stokes
python3 stokes_ss.py

# drivencavity steadystate
python3 drivencavity_steadystate.py

# drivencavity tdp pout vout distrcontrol
python3 drivencavity_tdp_pout_vout_distcontrol.py

# test units
python3 test_units.py
