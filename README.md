# 1D turbulent-driven wind outflow for Red Supergiant stars

This folder contains a simple hydrodynamic setup for modelling the effect of (constant) turbulent pressure on the wind outflow of a Red Supergiant star. It follows closely the steady-state model presented by [Kee et al. (2021), A&A 646](https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.180K/abstract). Some experimentation has been done with a spatially-decreasing turbulent pressure (as should be expected from turbulent decay), but this is not included in the present code. Namely, no viable wind solutions are obtained for somewhat realistic spatial variation of turbulent pressure.

It was developed for the purpose of learning more about winds from luminous, cool stars while working at the Institute of Astronomy of KU Leuven (2025).

## Setup

After cloning go into the local directory and fetch the AMRVAC makefile (assuming the `AMRVAC_DIR` variable is exported within the `.bashrc` on linux or `.bash_profile` on macos)
```
$ setup.pl -d=1
```
and do a (parallel) compilation and if the AMRVAC source was not compiled yet, this will be done now as well,
```
$ make -j
```

## How-to

### Additional user setup

I also use my own routine that contains all Fortran kind specifications. My external routines are not available yet (**TO DO**). But to ensure correct compilation remove the following declaration in the module header
```
use mod_kind_parameter, only: dp
```
and add instead
```
! Double precision kind
integer, parameter :: dp = kind(0.d0)
```
After this rename or delete the `local.make` file to not include it in the main compilation.

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the parameter file can be looked up on the AMRVAC webpage. The parameter file provided is 'basic' that works for the problem. Additionally, a `star_list` is specified in the `usr.par` file containing variables specific for this user problem. A `convert.par` file is included to post-process the binary output files to simple ASCII files.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar_sol    | stellar mass (solar units)                                    |
| rstar_sol    | stellar radius (solar units)                                  |
| twind_cgs    | wind temperature (Kelvin)                                     |
| vturb_cgs    | (constant) turbulent speed (cm/s) throughout wind             |
| kappa_ross_cgs | Rosseland mean opacity (cm^2/g)                             |

## Notice

Tested with AMRVAC version 3.2 and 3.1.

## Known issues

None
