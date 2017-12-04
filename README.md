## Example Setups of Incompressible Navier--Stokes Equations with Control and Observation as Linear-quadratic State Space Systems

This repo contains the data and example code to simulate incompressible flows with actuation and observation only using tools from a linear algebra package. The main features are

 * the (quadratic) nonlinearity is preassembled as a tensor
 * the Reynolds number dependent terms are separated and, thus, can be scaled to realize any Reynolds number
 * boundary control is realized through a Robin-type relaxation

### Tell me more

Please see the techreport [arxiv:1707.08711](https://arxiv.org/abs/1707.08711).

### Run the code

To get started, replicate the example setups.

#### Python

 1. Go to the `python` directory and create the folders needed for storing the computed data.
```
cd python
mkdir results
mkdir data
mkdir tikz
```
 2. Launch the provided example problem scripts. E.g.
```
python drivencavity_steadystate.py
```
 3. Check the output. In this case through
```
paraview results/v__drivencavity_stst.vtu
paraview results/p__drivencavity_stst.vtu
```

#### Octave/Matlab
 1. Go to the `matlab` directory and create the folders needed for storing the computed data.
```
cd matlab
mkdir results
```
 2. Launch `octave` or `matlab`, append the directories of the helper files to the runtime path, and run the provided example problem scripts. E.g.
```
octave
>> add_to_path
>> cylinderwake_steadystate
```
 3. Check the output. In this case -- with the default parameters -- through
```
paraview results/v__cylinderwake_stst_Re40_NV5812.vtu
paraview results/p__cylinderwake_stst_Re40_NV5812.vtu
```

### Example Outputs

Please see the techreport [arxiv:1707.08711](https://arxiv.org/abs/1707.08711).

### Support
The purpose, the concepts and the theoretical background of the code are laid out in the techreport. Feel free to contact the authors to discuss any conceptual issues.

This repo is a snapshot of `gitlab` which is under continuous development. To file or fix a bug or to request or add new features, please refer to there. 

### License
MIT

### Cite it
To acknowledge the use of the code or the data, please refer to the techreport
```
Bibtex
```

