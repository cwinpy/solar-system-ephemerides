import datetime
import os
import sys
import argparse
from typing import Union

from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time, TimeDelta
from astropy import constants as const

from ._version import __version__
from .ephemeris import BodyEphemeris

try:
    __author__ = os.environ["USER"]
except KeyError:
    __author__ = "Unknown"

__date__ = datetime.datetime.now()


HEADER = """\
# Build information for {exec}
# Author: {author}
# solar_system_ephemerides version: {version}
# ephemeris creation date: {date}
#
# Ephemeris creation command:-
#       {command}
#
# The JPL {ephem} ephemeris {ephemurl} has been used.
#
# This file consists of a header line containing:
#       GPS time of the start, interval between entries (secs), no. of entries
# Each entry consists of:
#       GPS time                Pos. x (lt sec)         Pos. y (lt sec)
#       Pos. z (lt sec)         Vel. x (lt sec/sec)     Vel. y (lt sec/sec)
#       Vel. z (lt sec/sec)     Acc. x (lt sec/sec^2)   Acc. y (lt sec/sec^2)
#       Acc. z (lt sec/sec^2)
"""

# set locations of JPL ephemeris files for downloading
EPH_URLS = {
    "DE440": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    "DE438": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de438.bsp",
    "DE436": "ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de436.bsp",
    "DE435": "ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de435.bsp",
    "DE432S": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp",
    "DE430": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp",
    "DE421": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp",
    "DE414": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de414.bsp",
    "DE410": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de410.bsp",
    "DE406S": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de406s.bsp",
    "DE405": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp",
    "DE200": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de200.bsp",
}

# code description
DESCRIPTION = """Create an ephemeris file for a given solar system body containing positions and \
velocities at a requested set of times.
"""

# set of allowed solar system bodies
BODIES = [
    "sun",
    "mercury",
    "venus",
    "earth-moon-barycenter",
    "earth-moon-barycentre",  # allow UK spelling
    "earth",
    "moon",
    "mars",
    "jupiter",
    "saturn",
    "uranus",
    "neptune",
    "pluto",
]


def cli():  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e",
        "--ephemeris",
        "--jplde",
        "--ephem",
        dest="ephemeris",
        required=True,
        help="Set the ephemeris to use (e.g. DE405)",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        dest="output",
        required=True,
        help="Set the output file. If this ends in '.gz' it will be gzipped.",
    )
    parser.add_argument(
        "-g",
        "--gps-start",
        dest="gpsstart",
        type=int,
        default=None,
        help="Set the GPS time at which to start generating the ephemeris",
    )
    parser.add_argument(
        "-y",
        "--year-start",
        dest="yearstart",
        type=int,
        default=None,
        help="Set the (decimal) year at which to start generating the ephemeris",
    )
    parser.add_argument(
        "-n",
        "--num-years",
        dest="nyears",
        type=float,
        required=True,
        help="Set the number of years over which to generate the ephemeris",
    )
    parser.add_argument(
        "-i",
        "--interval",
        dest="interval",
        type=float,
        required=True,
        help="Set the time step (in hours) between successive output points",
    )
    parser.add_argument(
        "-t",
        "--target",
        "--body",
        dest="target",
        required=True,
        help="Set the solar system body to generate the ephemeris for",
    )

    args = parser.parse_args()

    _ = generate_ephemeris(
        body=args.target,
        jplde=args.ephemeris,
        nyears=args.nyears,
        interval=args.interval,
        gpsstart=args.gpsstart,
        yearstart=args.yearstart,
        output=args.output,
        cli=True,
        overwrite=True,
    )


def generate_ephemeris(
    body: str = None,
    jplde: str = None,
    nyears: Union[int, float] = None,
    interval: Union[int, float] = None,
    gpsstart: int = None,
    yearstart: int = None,
    output: str = None,
    overwrite: bool = False,
    cli: bool = False,
) -> BodyEphemeris:
    """
    Generate an ephemeris for a given solar system body. Return a BodyEphemeris
    object.

    Parameters
    ----------
    body: str
        The solar system body for which to generate the ephemeris, e.g.,
        "earth".
    jplde: str
        The JPL development ephemeris version to use, e.g., "DE405".
    nyears: int, float
        The number of years over which to generate the ephemeris.
    interval: int, float
        The time step (in hours) between successive ephemeris values.
    gpsstart: int
        The start time for the ephemeris in GPS seconds.
    yearstart: int
        The start time for the ephemeris as a calendar year value.
    output: str
        The path to a file into which to output the ephemeris. If this ends
        with ".gz" the file will be gzipped.
    overwrite: bool
        If outputting to a file, set whether to overwrite an existing file or
        not. Default is False.
    """

    if jplde is None:
        raise ValueError("JPL development version must be given.")

    # check ephemeris is either a filename, or in our current list
    if jplde.endswith(".bsp"):  # pragma: no cover
        ephemfile = jplde
    else:
        if jplde.upper() not in EPH_URLS.keys():
            raise ValueError(
                f"Ephemeris '{jplde}' is not allowed, use one of: {EPH_URLS.keys()}."
            )
        else:
            ephemfile = EPH_URLS[jplde.upper()]

    # check that the body is in our current list
    if body.lower() not in BODIES:
        raise ValueError(f"Target body '{body}' is not in the allowed list: {BODIES}.")
    else:
        body = (
            body.lower()
            if body.lower() != "earth-moon-barycentre"
            else "earth-moon-barycenter"
        )

    # set the ephemeris file
    solar_system_ephemeris.set(ephemfile)

    # set the start time
    if gpsstart is not None and yearstart is not None:
        raise ValueError("Specify either '--gps-start' or '--year-start', but not both")
    elif gpsstart is not None:
        try:
            starttime = Time(gpsstart, format="gps", scale="utc")
        except Exception as e:
            raise ValueError(f"Could not parse start GPS time {gpsstart}: {e}")
    else:
        try:
            starttime = Time(yearstart, format="decimalyear", scale="utc")
        except Exception as e:
            raise ValueError(f"Could not parse start year {yearstart}: {e}")

    # set the time step
    try:
        dt = TimeDelta(interval * 3600.0, format="sec")
    except Exception as e:
        raise ValueError(f"Could not parse time interval {interval}: {e}")

    # set the end time
    try:
        endtime = (
            Time(starttime.decimalyear + nyears, format="decimalyear", scale="utc") + dt
        )
    except Exception as e:
        raise ValueError(f"Could not parse total timespan {nyears}: {e}")

    pos = []
    vel = []
    acc = []

    # get positions, velocities and accelerations
    curtime = starttime
    while curtime < endtime:
        tpos, tvel = get_body_barycentric_posvel(body, curtime)

        # convert positions to light seconds
        pos.append(tpos.xyz.to("m") / const.c)

        # convert velocities to light seconds per second
        vel.append(tvel.xyz.to("m/s") / const.c)

        # calculate accelerations (using velocities +/- dt/2 around the current time)
        ctminus = curtime - dt / 2.0
        _, mvel = get_body_barycentric_posvel(body, ctminus)
        ctplus = curtime + dt / 2.0
        _, pvel = get_body_barycentric_posvel(body, ctplus)
        acc.append(((pvel.xyz.to("m/s") - mvel.xyz.to("m/s")) / const.c) / dt.to("s"))

        curtime += dt

    # store and return as BodyEphemeris object
    be = BodyEphemeris()
    be.pos = [tpos.value for tpos in pos]
    be.vel = [tvel.value for tvel in vel]
    be.acc = [tacc.value for tacc in acc]

    be.body = body
    be.t0 = starttime.gps
    be.timestep = dt.value
    be.nentries = len(pos)

    if output is not None:
        # set output header
        headdic = {}
        headdic["exec"] = os.path.basename(sys.argv[0]) if cli else "N/A"
        headdic["author"] = __author__
        headdic["date"] = __date__
        headdic["version"] = __version__
        headdic["command"] = (
            " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]) if cli else "N/A"
        )
        headdic["ephem"] = jplde.upper()
        headdic["ephemurl"] = ephemfile

        be.write(output, header=HEADER.format(**headdic), overwrite=overwrite)

    return be
