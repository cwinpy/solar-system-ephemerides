# Solar System ephemerides

A package holding solar system ephemeris files, storing positions, velocities and accelerations of
the Earth and Sun for a range of [JPL development
ephemeris](https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris) versions.
These can be used, for example, for calculating Doppler modulations and relativistic corrections for
continuous gravitational-wave signals. The package contains a command line script and Python API for
generating these files, which makes use of the [Astropy](https://www.astropy.org/) and
[jplephem](https://pypi.org/project/jplephem/) packages. The script can more generally be used to
create ephemeris files for any major planet within the solar system.

Along with ephemerides for the Earth and Sun, the package also contains two files providing time
corrections between [Terrestrial Time](https://en.wikipedia.org/wiki/Terrestrial_Time) (TT) and
[Barycentric Coordinate Time](https://en.wikipedia.org/wiki/Barycentric_Coordinate_Time) (TCB) or
[Barycentric Dynamical Time](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time) (TDB). These
were created using a script within
[LALSuite](https://lscsoft.docs.ligo.org/lalsuite/lalsuite/index.html), called
`lalpulsar_create_time_correction_ephemeris`, which uses the
[`TIMEEPH_short.te405`](https://bitbucket.org/psrsoft/tempo2/src/master/T2runtime/ephemeris/TIMEEPH_short.te405)
(for TCB conversions)
[`TDB.1950.2050`](https://bitbucket.org/psrsoft/tempo2/src/master/T2runtime/ephemeris/TDB.1950.2050)
(for TDB conversions) within [Tempo2](https://bitbucket.org/psrsoft/tempo2).

## Installation

This package can be installed, along with its requirements, from
[PyPI](https://pypi.org/project/solar-system-ephemerides) via `pip` with:

```bash
pip install solar-system-ephemerides
```

or, within a conda environment, from
[conda-forge](https://anaconda.org/conda-forge/solar_system_ephemerides) with:

```bash
conda install -c conda-forge solar_system_ephemerides
```

It can be installed from source with:

```bash
git clone git@github.com:cwinpy/solar-system-ephemerides.git
cd solar-system-ephemerides
pip install .
```

## Usage

Once installed, to get the path to an ephemeris file, e.g., the DE405 file for Earth, within Python,
you can do:

```python
from solar_system_ephemerides import body_ephemeris_path

path = body_ephemeris_path("earth", DE405)
```

Or, the paths can be output using the `ephemeris_path` command line script, e.g., with:

```bash
ephemeris_path --body earth --ephem DE405
```

The full usage of `ephemeris_path` is shown with the `--help` flag:

```bash
$ ephemeris_path --help
usage: ephemeris_path [-h] [--body [BODY [BODY ...]]]
                      [--ephem [EPHEM [EPHEM ...]]]
                      [--units [UNITS [UNITS ...]]] [--return-dir]
                      [--rel-path REL_PATH]

Return the path to an ephemeris file or directory containing ephemeris files

optional arguments:
  -h, --help            show this help message and exit
  --body [BODY [BODY ...]], -b [BODY [BODY ...]]
                        The solar system body[ies] to return the path for.
                        These must be in "earth" or "sun".
  --ephem [EPHEM [EPHEM ...]], -e [EPHEM [EPHEM ...]]
                        The JPL development ephemeris version(s) to use. These
                        must be in "DE200", "DE405", "DE421", "DE430", "DE435"
                        or "DE436".
  --units [UNITS [UNITS ...]], -u [UNITS [UNITS ...]]
                        The time system units to return the path for. This
                        must be "TCB" or "TDB".
  --return-dir, -d      Return the ephemeris directory for a given body/time
                        system unit rather than the file path.
  --rel-path REL_PATH, -r REL_PATH
                        Return paths relative to this given path.
```

which also shows the JPL development ephemeris versions and solar system bodies for which there are
ephemeris files available within the package.

## Ephemeris generation

The package come with a script called `create_solar_system_ephemeris` that can be used to generate
new ephemeris files.

> Note: this script is based on the
> [LALSuite](https://lscsoft.docs.ligo.org/lalsuite/lalsuite/index.html) script,
> `lalpulsar_create_solar_system_ephemeris_python`, and is designed to supersede it and the
> equivalent `C` script `lalpulsar_create_solar_system_ephemeris`.

For example, to create an ephemeris file for the Sun using the JPL DE421
spanning from 2015 to 2035, with 20 hours between each value, one would do:

```bash
create_solar_system_ephemeris --target SUN --year-start 2015 --interval 20 --num-years 20 --ephemeris DE421 --output-file sun15-35-DE421.dat.gz
```

Full usage of the script can be shown with the `--help` flag:

```bash
$ create_solar_system_ephemeris --help
usage: create_solar_system_ephemeris [-h] -e EPHEMERIS -o OUTPUT [-g GPSSTART]
                                     [-y YEARSTART] -n NYEARS -i INTERVAL -t
                                     TARGET

optional arguments:
  -h, --help            show this help message and exit
  -e EPHEMERIS, --ephemeris EPHEMERIS, --jplde EPHEMERIS, --ephem EPHEMERIS
                        Set the ephemeris to use (e.g. DE405)
  -o OUTPUT, --output-file OUTPUT
                        Set the output file. If this ends in '.gz' it will be
                        gzipped.
  -g GPSSTART, --gps-start GPSSTART
                        Set the GPS time at which to start generating the
                        ephemeris
  -y YEARSTART, --year-start YEARSTART
                        Set the (decimal) year at which to start generating
                        the ephemeris
  -n NYEARS, --num-years NYEARS
                        Set the number of years over which to generate the
                        ephemeris
  -i INTERVAL, --interval INTERVAL
                        Set the time step (in hours) between successive output
                        points
  -t TARGET, --target TARGET, --body TARGET
                        Set the solar system body to generate the ephemeris
                        for
```

[![PyPI version](https://badge.fury.io/py/solar_system_ephemerides.svg)](https://badge.fury.io/py/solar_system_ephemerides) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/solar_system_ephemerides.svg)](https://anaconda.org/conda-forge/solar_system_ephemerides)
