# OpenSim plugin for calculating the feasible muscle forces

## Description

This project contains an OpenSim v4.0 plugin for calculating the feasible muscle
forces that satisfy the movement and physiological muscle constraints. The
plugin also contains two executables, namely one for calculating the feasible
muscle forces and one for calculating the feasible joint reaction loads.

- *data* contains OpenSim .osim model related files
- *scripts* python scripts for plotting results

## Dependencies

This project has been tested on Linux. On Windows there is an issue with the gmp
library, which is used for the vertex enumeration (lrs).

- OpenSim v4.0 C++ API
- OpenSim v4.0 Python 3.7 bindings
- gmp library (https://gmplib.org/)

## Acknowledgments

[1] D. Stanev and K. Moustakas, Modeling Musculoskeletal Kinematic and Dynamic
Redundancy Using Null Space Projection, PLoS ONE, 14(1): e0209171, Jan. 2019,
https://doi.org/10.1371/journal.pone.0209171

[2] D. Stanev and K. Moustakas, Stiffness Modulation of Redundant
Musculoskeletal Systems, Journal of Biomechanics, vol. 85, pp. 101-107,
Mar. 2019, https://doi.org/10.1016/j.jbiomech.2019.01.017

[3] SimTK project: https://simtk.org/projects/redundancy

[4] SimTK project: https://simtk.org/projects/stiffness

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img
alt="Creative Commons License" style="border-width:0"
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is
licensed under a <a rel="license"
href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution
4.0 International License</a>.
