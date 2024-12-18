# A Pseudo-3D-Stokes Module for Varying Apertures
[![Identifier](https://img.shields.io/badge/doi-10.18419%2Fdarus--4313-d45815.svg)](https://doi.org/10.18419/darus-4313)
[![Identifier](https://img.shields.io/badge/Publication-blue)]([http://ssrn.com/abstract=4927521](https://doi.org/10.1016/j.advwatres.2024.104860))

**Pseudo-3D-Stokes Module for Varying Apertures** is a [DuMu<sup>x</sup>](https://dumux.org/) module developed at research institutions. [DuMu<sup>x</sup>](https://dumux.org/) is a simulation framework focusing on Finite Volume discretization methods, model coupling for multi-physics applications, and flow and transport applications in porous media.

This module aims to assist researchers in planning, improving, or interpreting microfluidic experiments through numerical simulations based on the Stokes equations in an easy and intuitive way. It uses `.pgm` files as input to create numerical grids. These `.pgm` files should include 8 bit grayscale values referring to the relative height of a microfluidic cell, which can be created from microscopy images of a microfluidic experiment using suitable image processing procedures.

Based on the `.pgm` files, the `python` module `localdrag` ( [git](https://git.iws.uni-stuttgart.de/krachdd/localdrag), [DaRUS](https://doi.org/10.18419/darus-4313) ) should be used to create the suitable drag prefactor fields `lambda1` and `lambda2`. For further details, refer to our publication [**A Novel Geometry-Informed Drag Term Formulation for Pseudo-3D Stokes Simulations with Varying Apertures**](http://ssrn.com/abstract=4927521).

## How to Install

The **Pseudo-3D-Stokes Module for Varying Apertures** is based on DuMu<sup>x</sup>. Therefore, DuMu<sup>x</sup> needs to be installed on a Linux system or within a Docker container. Instructions for installing DuMu<sup>x</sup> can be found [here](https://dumux.org/docs/doxygen/master/installation.html). Please clone the `master` branch.

In addition to DuMu<sup>x</sup> and the Dune core modules, the Dune module `subgrid` needs to be installed using the `installexternal.py` script provided within the DuMu<sup>x</sup> installation.

Next, clone the **Pseudo-3D-Stokes Module for Varying Apertures** into the DuMu<sup>x</sup> installation directory with the following command:
```bash
# with https 
git clone --recurse-submodules https://git.iws.uni-stuttgart.de/krachdd/pseudo3D_stokes.git
# with ssh key
git clone --recurse-submodules git@git.iws.uni-stuttgart.de:krachdd/pseudo3D_stokes.git
```
`localdrag` is included via die `--recurse-submodules` flag as a submodule and cloned into the folder `preprocessing`. Please check for python requirements in `cat pseudo3D_stokes/preprocessing/localdrag/requirements.txt`. 
Checkout main branch and add the `localdrag` module's directory to the path directory list using `PYTHONPATH`
```bash
cd pseudo3D_stokes 
# get the current main branch of the submodule
git submodule foreach --recursive git checkout main
# export the absolute path 
export PYTHONPATH=$PYTHONPATH:$(pwd)/preprocessing/localdrag/
# echo the path - should not be empty!
echo $PYTHONPATH
cd .. 
```
Subsequently the dunecontrol command needs to be executed:
```bash
# Reconfiguration and build of all modules
./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all
```
The ctests can be executed to make sure everything has been installed properly:
```bash
# go to build directory
cd pseudo3D_stokes/build-cmake/appl/
# make executable
make pseudo3D_stokes
# run all 8 tests
ctest all
```
## How to run simulations
To start, it is recommended to run the prepared test case **test_singlePrecipitate_2d_total**. This folder contains `.pgm` files representing a domain with a semispherical precipitate at the center. The varying grayscale values in these files represent the apertures of the microfluidic cell (0: no aperture; 255: maximum aperture). Additionally, for each `.pgm` file, there are two `.txt` files with the prefixes `lambda1` and `lambda2`, which include the local drag terms for each grid point. These drag terms are calculated using the preprocessing `python` module [`localdrag`](https://git.iws.uni-stuttgart.de/krachdd/localdrag), based on the `.pgm` files. Further details about the preprocessing step can be found in the [`localdrag` module documentation](https://git.iws.uni-stuttgart.de/krachdd/localdrag/-/blob/faa44b05c1296d94a6f2fedc9752105ab20ef511/README.md).
In order to run the test case, you need to first change to the correct directory, and build the executable - analogously to the last step in the *How to install*: 
```bash
cd pseudo3D_stokes/build-cmake/appl/
make pseudo3D_stokes
```
In this folder, the `python` script **_runStokesGeneric.py_** is located. This script executes multiple simulations based on a defined _directory_, _height of the microfluidic channel_, and _voxel size_. You need to provide these input parameters when running the script. To learn more about the possible input parameters and how to use them, run the following command:

```bash
python3 runStokesGeneric.py --help
```
This command will display a help message with detailed information on the available options and usage instructions.

We recommend using the shell script [run_simulations.sh](appl/run_simulations.sh) contains a command line that executes the `python` script [runStokesGeneric.py](appl/runStokesGeneric.py) with predefined input parameters. To run the simulations for the test case with these predefined parameters, use the following command:
```bash
sh run_simulations.sh
```
To view the simulation results in ParaView, use:
```bash
paraview test_singlePrecipitate_2d_total/*.pvd&
```

If you want to run simulations with your own geometries, follow these steps:

1.  **Create a new directory** containing your `.pgm` files along with the corresponding lambda files. These lambda files can be generated using the preprocessing `python` module [`localdrag`](https://git.iws.uni-stuttgart.de/krachdd/localdrag).
    
2.  **Run the simulations** by adjusting the input parameters in [runStokesGeneric.py](appl/runStokesGeneric.py) or better using [run_simulations.sh](appl/run_simulations.sh) after configuring it for your new setup.


## How to cite

If you are using the **Pseudo-3D-Stokes module for varying apertures** module in scientific publications and in
the academic context, please cite our publication:
**A novel geometry-informed drag term formulation for pseudo-3D Stokes simulations with varying apertures.**

```bib
@article{Krach2025,
    title = {A novel geometry-informed drag term formulation for pseudo-3D Stokes simulations with varying apertures},
    journal = {Advances in Water Resources},
    volume = {195},
    year = {2025},
    doi = {https://doi.org/10.1016/j.advwatres.2024.104860},
    author = {David Krach and Felix Weinhardt and Mingfeng Wang and Martin Schneider and Holger Class and Holger Steeb},
    keywords = {Porous media, Stokes flow, Biomineralization, Microfluidics, Image-based simulations, Computational efficiency versus accuracy}
}
```

## Acknowledgements
Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy (Project number 390740016 - EXC 2075 and the Collaborative Research Center 1313 (project number 327154368 - SFB1313). We acknowledge the support by the Stuttgart Center for Simulation Science (SimTech).

## Developer
- [Felix Weinhardt](https://www.mib.uni-stuttgart.de/de/institut/team/Weinhardt-00003/) E-mail: [felix.weinhardt@mib.uni-stuttgart.de](mailto:felix.weinhardt@mib.uni-stuttgart.de)
- [David Krach](https://www.mib.uni-stuttgart.de/institute/team/Krach/) E-mail: [david.krach@mib.uni-stuttgart.de](mailto:david.krach@mib.uni-stuttgart.de)

## Contact
- [Software Support Institute of Applied Mechanics](mailto:software@mib.uni-stuttgart.de)
- [Data Support Institute of Applied Mechanics](mailto:data@mib.uni-stuttgart.de)

