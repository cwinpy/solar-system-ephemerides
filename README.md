# Solar System ephemerides

A package holding solar system ephemeris files, storing positions and velocities of the Earth and
Sun for a range of [JPL development
ephemeris](https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory_Development_Ephemeris) versions.
These can be used, for example, for calculating Doppler modulations for continuous
gravitational-wave signals.

## Installation

This package can be installed, along with its requirements, from PyPI via `pip` with:

```bash
pip install solar_system_ephemerides
```

It can be installed from source with:

```bash
git clone git@git.ligo.org:CW/ephemerides/solar-system-ephemerides.git
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

## Ephemeris generation

The package come with a script called `create_solar_system_ephemeris` that can be used to generate
new ephemeris files. For example, to create an ephemeris file for the Sun using the JPL DE421
spanning from 2015 to 2035, with 20 hours between each value, one would do:

```bash
create_solar_system_ephemeris --target SUN --year-start 2015 --interval 20 --num-years 20 --ephemeris DE421 --output-file sun15-35-DE421.dat.gz
```

[![PyPI version](https://badge.fury.io/py/solar_system_ephemerides.svg)](https://badge.fury.io/py/solar_system_ephemerides)
