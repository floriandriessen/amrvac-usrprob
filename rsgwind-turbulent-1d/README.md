# 1D turbulent-driven wind outflow for Red Supergiant stars

This folder contains a simple hydrodynamic setup for modelling the effect of (constant) turbulent pressure on the wind outflow of a Red Supergiant star. It follows closely the steady-state model presented by [Kee et al. (2021), A&A 646](https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.180K/abstract). Some experimentation has been done with a spatially-decreasing turbulent pressure (as should be expected from turbulent decay), but this is not included in the present code. Namely, no viable wind solutions are obtained for somewhat realistic spatial variation of turbulent pressure (i.e. exponential decay over a representative scale height).

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
A `local.make` file is included that ensures external user modules from the main directory `share/` folder are visible during compilation. The main makefile calls `local.make` and will use its information for the compilation of `mod_usr.t`--and thus the final executable.

## How-to

### Performing a simulation

Simulation follows the "typical" model of Kee et al. (2021, sec. 3.1). For example, a run on a single core
```
$ ./amrvac -i usr.par
```
and once finished to convert a range of output
```
aiconvert convert.par N0 N1
```
where $N0$, $N1$ are an integer for the desired initial and final snapshot. *If changes occur in the namelist(s) meshlist and/or usr_list, this corresponding change also needs to be made in the convert.par file.*

**These models should only be run on a single core.** This is because the ray tracing done to compute the optical depth, which in turn is used to compute the total thermal pressure (gas + turbulent), is outside-in oriented. In parallel runs there is no guarantee that the block ordering is correctly aligned due to the Morton curve (even more a problem in multi-D) and so the optical depth might be computed erroneously.

Still given the restriction the simulations are rather fast (few minutes depending on grid setup) so there are no options included in the .par file to do restarts or resumes. Depending on the stellar parameters the runtime can be best adjusted based on the `unit_time` variable that is printed to the screen at the start of the simulation.

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the parameter file can be looked up on the AMRVAC webpage. The parameter file provided is basic that works for the problem. Additionally, a `star_list` is specified in the `usr.par` file containing variables specific for this user problem. A `convert.par` file is included to post-process the binary output files to simple ASCII files.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar_sol    | stellar mass (solar units)                                    |
| rstar_sol    | stellar radius (solar units)                                  |
| twind_cgs    | wind temperature (Kelvin)                                     |
| vturb_cgs    | (constant) turbulent speed (cm/s) throughout wind             |
| kappa_ross_cgs | Rosseland mean opacity (cm^2/g)                             |

## Notice

Tested with AMRVAC version 3.2 and 3.1.

*Backward compatibility of version 3.x with older versions (especially < v3.x) is not the case due to internal AMRVAC code changes. Within a given main version number backward compatibility **might** be guaranteed--always consult the AMRVAC release log if notable changes occurred.*

## Known issues

- The first physical cell seems to always attain a slightly negative velocity--no matter what BC is taken, initial condition, numerical scheme combination, etc. This does not happen in the steady-state model of Kee et al. Nonetheless, it has no effect on the final dynamical wind solution and mass-loss rate.
