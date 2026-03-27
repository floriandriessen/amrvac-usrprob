
# 1-D dusty AGB pulsation-driven wind

A simple model to study the dusty pulsation-driven winds of AGB stars as first pioneered by [Bowen (1988)](https://ui.adsabs.harvard.edu/abs/1988ApJ...329..299B/abstract). The dust implementation follows Bowen's phenomenological model and is further simplified by assuming a power-law radiative equilibrium temperature structure. This is contrary to Bowen where the radiative equilibirum temperature is computed from the local spherically-modified optical depth. The difference appears, however, of minor importance in model results.

## Setup

After cloning go into the local directory and fetch the makefile (assuming you have exported the `AMRVAC_DIR` variable in the `.bashrc` on Linux or `.bash_profile` on macOS)
```
$ setup.pl -d=1 -v=1
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```
A `local.make` file is included that ensures external user modules from the main directory `share/` folder are visible during compilation. The main makefile calls `local.make` and will use its information for the compilation of `mod_usr.t`--and thus the final executable.

## How-to

### Performing a simulation

Simulation follows the default dust model setup of Bowen (1988). For example, a run on 4 cores
```
$ mpiexec -np 4 ./amrvac -i usr.par
```
and once finished to convert a range of output
```
aiconvert convert.par N0 N1
```
where $N0$, $N1$ are an integer for the desired initial and final snapshot. *If changes occur in the namelist(s) meshlist and/or usr_list, this corresponding change also needs to be made in the convert.par file.*

Given that the simulations are rather fast (few seconds to few minutes depending on grid setup) there are no options included in the .par file to do restarts or resumes. Depending on the stellar parameters the runtime can be best adjusted based on the `unit_time` variable that is printed to the screen at the start of the simulation.

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are basic settings that work for the problem.

Additionally, a `usr_list` is specified in the .par file containing variables specific for our problem. The ones relevant for computations are converted in the code to dimensionless quantities within the `initglobaldata_usr` subroutine.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar_sol    | stellar mass (solar units)                                    |
| v_amp_cgs    | velocity amplitude of pulsation (cm/s)                        |
| puls_period_day | pulsation period (days)                                    |
| teff_cgs     | stellar effective temperature (K)                             |
| rhostar_cgs  | stellar boundary density (g/cm^3)                             |
| beta         | power-law index for radiative equilibrium temperature stratification |
| tcond1__cgs  | condensation temperature of dust (K)                          |
| grad_to_ggrav | ratio of radiation to gravity force                          |
| use_lya_cooling | activate Lyman-alpha cooling                               |

The stellar radius is computed from the Period-Mass relation of Ostlie & Cox applicable to Mira stars. The stellar luminosity then follows directly with the given input effective temperature.

## Notice

Tested with AMRVAC version 3.3, 3.2. 

*Backward compatibility of version 3.x with older versions (especially < v3.x) is not the case due to internal AMRVAC code changes. Within a given main version number backward compatibility **might** be guaranteed--always consult the AMRVAC release log if notable changes occurred.*

## Known issues

None.
