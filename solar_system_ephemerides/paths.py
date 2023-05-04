import gzip
import sys
from argparse import ArgumentParser
from pathlib import Path

import pkg_resources


EPHEMERIS_BASE_PATH = pkg_resources.resource_filename(
    "solar_system_ephemerides", "ephemerides/"
)

# dictionary of current bodies stored within the package (with aliases)
BODIES = {
    "earth": ["earth", "terra", "terre", "erde", "tierra", "ziemia", "aarde"],
    "sun": ["sun", "sol", "sonne", "zon", "sole", "słońce"],
}

# list of valid JPL DE values stored within the package
JPLDE = [
    "DE200",
    "DE405",
    "DE421",
    "DE430",
    "DE435",
    "DE436",
]

# list of valid timespans for ephemerides
TIMESPANS = [
    "00-40",
]

# dictionary of time correction files prefixes
TIMECORRS = {
    "TDB": "tdb_2000-2040.dat.gz",
    "TCB": "te405_2000-2040.dat.gz",
}


class EphemerisPath:
    def __init__(self, body: str, jplde: str, timespan: str = "00-40"):
        self.body = body
        self.jplde = jplde
        self.timespan = timespan

    @property
    def body(self):
        return self._body

    @body.setter
    def body(self, body):
        if not isinstance(body, str):
            raise TypeError("body must be a string")

        if body.lower() in BODIES:
            self._body = body.lower()
        else:
            for b in BODIES:
                if body.lower() in BODIES[b]:
                    self._body = body.lower()
                    break
            else:
                raise ValueError(
                    f"body '{body}' is not recognised. It must be in {list(BODIES.keys())}."
                )

    @property
    def jplde(self):
        return self._jplde

    @jplde.setter
    def jplde(self, jplde):
        destr = str(jplde).upper()

        # check if string starts with "DE", if not prepend "DE"
        if not destr.startswith("DE"):
            destr = "DE" + destr

        if destr in JPLDE:
            self._jplde = destr
        else:
            raise ValueError(f"JPL DE must be in {JPLDE}")

    @property
    def timespan(self):
        return self._timespan

    @timespan.setter
    def timespan(self, timespan):
        if isinstance(timespan, (tuple, list)):
            if len(timespan) == 2:
                syear = str(timespan[0])[-2:]  # truncate, e.g, 2000 to just 00
                eyear = str(timespan[1])[-2:]
                tstr = f"{syear}-{eyear}"
            else:
                raise ValueError("timespan must contain start year and end year")
        elif isinstance(timespan, str):
            tstr = timespan
        else:
            raise TypeError("timespan must be a string or pair of integer")

        if tstr in TIMESPANS:
            self._timespan = tstr
        else:
            raise ValueError(f"time span must be in {TIMESPANS}")

    @property
    def path(self):
        """
        Construct the path. Return a pathlib.Path object.
        """

        path = Path(EPHEMERIS_BASE_PATH)
        path = path / self.body / f"{self.body}{self.timespan}-{self.jplde}.dat.gz"

        return path

    def __call__(self):
        return self.path

    def __str__(self):
        return str(self.path)

    @property
    def contents(self):
        """
        Return the contents of the file in a string.
        """

        with gzip.open(self.path, "r") as fp:
            contents = fp.read()

        return contents


class TimeEphemerisPath:
    def __init__(self, units: str):
        self.units = units

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        if not isinstance(units, str):
            raise TypeError("units must be a string")

        if units.upper() in TIMECORRS:
            self._units = units.upper()
        else:
            raise ValueError(
                f"units '{units}' is not recognised. It must be in {list(TIMECORRS.keys())}."
            )

    @property
    def path(self):
        """
        Construct the path. Return a pathlib.Path object.
        """

        path = Path(EPHEMERIS_BASE_PATH)
        path = path / "time" / f"{TIMECORRS[self.units]}"

        return path

    def __call__(self):
        return self.path

    def __str__(self):
        return str(self.path)

    @property
    def contents(self):
        """
        Return the contents of the file in a string.
        """

        with gzip.open(self.path, "r") as fp:
            contents = fp.read()

        return contents


def ephemeris_path(
    body: str, jplde: str = "DE405", timespan: str = "00-40", string: bool = False
) -> str:
    """
    Function for returning an ephemeris file path.

    Parameters
    ----------
    body: str
        The solar system body to get the ephemeris for, e.g., "earth".
    jplde: str
        The JPL development ephemeris version, e.g., "DE405".
    timespan: str
        The ephemeris timespan. Currently only 00-40 is included, so this is the default value.
    string: bool
        If True, return the path as a string rather than a pathlib.Path object.
    """

    path = EphemerisPath(body=body, jplde=jplde, timespan=timespan)

    return str(path) if string else path()


def time_ephemeris_path(units: str, string: bool = False):
    """
    Function for returning a time correction ephemeris file path.

    Parameters
    ----------
    units: str
        The time correction file units, either "TDB" or "TCB"
    string: bool
        If True, return the path as a string rather than a pathlib.Path object.
    """

    path = TimeEphemerisPath(units=units)

    return str(path) if string else path()


def cli():
    """
    Entry point for command line interface for returning paths.
    """

    parser = ArgumentParser(
        description=(
            "Return the path to an ephemeris file or directory containing "
            "ephemeris files"
        ),
    )

    bodies = ", ".join([f'"{body}"' for body in list(BODIES.keys())[:-1]])
    bodies += f' or "{list(BODIES.keys())[-1]}"'

    ephems = ", ".join([f'"{ephem}"' for ephem in JPLDE[:-1]])
    ephems += f' or "{JPLDE[-1]}"'

    parser.add_argument(
        "--body",
        "-b",
        default=None,
        nargs="*",
        help=f"The solar system body[ies] to return the path for. These must be in {bodies}.",
    )
    parser.add_argument(
        "--ephem",
        "-e",
        default=None,
        nargs="*",
        help=f"The JPL development ephemeris version(s) to use. These must be in {ephems}.",
    )
    parser.add_argument(
        "--units",
        "-u",
        default=None,
        nargs="*",
        help='The time system units to return the path for. This must be "TCB" or "TDB".',
    )
    parser.add_argument(
        "--return-dir",
        "-d",
        action="store_true",
        default=False,
        help="Return the ephemeris directory for a given body/time system unit rather than the file path.",
    )

    args = parser.parse_args()

    paths = []

    if args.body is not None:
        for body in args.body:
            ephemp = args.ephem if args.ephem is not None else ([JPLDE[0]] if args.return_dir else [])
            for ephem in ephemp:
                try:
                    paths.append(ephemeris_path(body=body, jplde=ephem))
                except (TypeError, ValueError):
                    sys.exit(1)

    if args.units is not None:
        for unit in args.units:
            try:
                paths.append(time_ephemeris_path(units=unit))
            except (TypeError, ValueError):
                sys.exit(1)

    if len(paths):
        if args.return_dir:
            print(*list(set([str(path.parent) for path in paths])), sep=":")
        else:
            print(*[str(path) for path in paths], sep=":")

        sys.exit(0)
    else:
        sys.exit(1)