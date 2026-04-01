
# 2.5D Alfven-wave-driven wind

This user code is an exploration to better understand Alfven-wave-driven winds from the Sun and other stars. It aims to reproduce some of the results presented in the papers of [Ofman & Davila (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJ...476..357O/abstract), applied to the Sun, and [Airapetian et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...528..965A/abstract) who applied the Ofman & Davila model to Red Supergiants.

It was developed for the purpose of learning more about winds from luminous, cool stars while working at the Institute of Astronomy of KU Leuven (2025).

## Setup

After cloning go into the local directory and fetch the AMRVAC makefile (assuming the `AMRVAC_DIR` variable is exported within the `.bashrc` on linux or `.bash_profile` on macos)
```
$ setup.pl -d=2 -v=3
```
and do a (parallel) compilation and if the AMRVAC source was not compiled yet, this will be done now as well,
```
$ make -j
```
A `local.make` file is included that ensures external user modules from the main directory `share/` folder are visible during compilation. The main makefile calls `local.make` and will use its information for the compilation of `mod_usr.t`--and thus the final executable.

## How-to

### Performing a simulation

The model setup of Ofman & Davila or Airapetian et al. can be chosen as reference by the name of the parameter file. For example,
```
$ mpiexec -np 4 ./amrvac -i ofman-davila1997.par
```
and once finished to convert a range of output
```
aiconvert convert-ofman-davila1997.par N0 N1
```
where $N0$, $N1$ are an integer for the desired initial and final snapshot. *If changes occur in the namelist(s) meshlist and/or usr_list, this corresponding change also needs to be made in the convert.par file.*

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the parameter file can be looked up on the AMRVAC webpage. The parameter file provided is basic that works for the problem. Additionally, a `star_list` is specified in the `usr.par` file containing variables specific for this user problem. A `convert.par` file is included to post-process the binary output files to simple ASCII files.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar_sol    | stellar mass (solar units)                                    |
| rstar_sol    | stellar radius (solar units)                                  |
| twind_cgs    | wind temperature (Kelvin)                                     |
| b0_cgs       | stellar magnetic field strength (Gauss)                       |
| rho0_cgs     | stellar surface mass density (g/cm^3)                         |
| vdrive_cgs   | perturbation amplitude of Alfven wave (cm/s)                  |

## Notice

Tested with AMRVAC version 3.3, 3.2.

*Backward compatibility of version 3.x with older versions (especially < v3.x) is not the case due to internal AMRVAC code changes. Within a given main version number backward compatibility **might** be guaranteed--always consult the AMRVAC release log if notable changes occurred.*

## Known issues

- The current user problem eventually crashes due to some infall of matter coming from the outer radial grid. This infall happens gradually and the hydrodynamic timestep slowly decreases up to a point it suddenly explodes to a very high timestep. The papers do not suffer from this and the reason of occurrence here remains unclear. The Ofman & Davila setup works better, that is more stable, than the Airapetian setup.
 

