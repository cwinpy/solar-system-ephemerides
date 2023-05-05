import numpy as np
from gzip import open as gzopen

from .paths import body_ephemeris_path, time_ephemeris_path


def decode(x):
    # decode bytes
    try:
        return x.decode()
    except AttributeError:
        return x


class BodyEphemeris:
    def __init__(self, body=None, jplde=None, timespan=None):
        self.body = body
        self.jplde = jplde
        self.timespan = timespan

        # try getting ephemeris from a file
        self.get_ephemeris()

    def get_ephemeris(self, body=None, jplde=None, timespan=None):
        if body is not None:
            # set body to given body
            self.body = body

        if jplde is not None:
            self.jplde = jplde

        try:
            ts = (
                timespan
                if timespan is not None
                else (self.timespan if self.timespan is not None else "00-40")
            )
            path = body_ephemeris_path(body=self.body, jplde=self.jplde, timespan=ts)
            self.timespan = ts
        except (ValueError, TypeError):
            # exit without creating an ephemeris
            return None

        # read in file
        try:
            fopen = gzopen if str(path).endswith(".gz") else open
            with fopen(path, "r") as fp:
                # data lines with contain 10 values
                data = np.genfromtxt(
                    [
                        decode(line)
                        for line in fp.readlines()
                        if len(decode(line).split()) == 10
                    ],
                    comments="#",
                )
        except Exception as e:
            raise IOError(f"Could not read ephemeris file: {e}")

        # set properties
        self.t0 = data[0, 0]  # start time of ephemeris (GPS time)
        self.timestep = data[1, 0] - self.t0  # time step (s)
        self.nentries = data.shape[0]  # number of entries

        self.pos = data[:, 1:4]
        self.vel = data[:, 4:7]
        self.acc = data[:, 7:10]

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, pos):
        self._pos = np.asarray(pos, dtype=float)

        if self.pos.shape[-1] != 3:
            raise ValueError("Position must contain x, y and z components")

    @property
    def vel(self):
        return self._vel

    @vel.setter
    def vel(self, vel):
        self._vel = np.asarray(vel, dtype=float)

        if self.vel.shape[-1] != 3:
            raise ValueError("Velocity must contain x, y and z components")

    @property
    def acc(self):
        return self._acc

    @acc.setter
    def acc(self, acc):
        self._acc = np.asarray(acc, dtype=float)

        if self.acc.shape[-1] != 3:
            raise ValueError("Acceleration must contain x, y and z components")

    @property
    def times(self):
        if (
            hasattr(self, "t0")
            and hasattr(self, "timestep")
            and hasattr(self, "nentries")
        ):
            return np.arange(
                self.t0, self.t0 + self.timestep * self.nentries, self.nentries
            )
        else:
            return None

    def __len__(self):
        if hasattr(self, "nentries"):
            return self.nentries
        else:
            return 0


def lal_ephemeris_data(jplde: str = "DE405"):
    """
    Return a LAL EphemerisData object containing the ephemerides for the Earth
    and Sun.

    Parameters
    ----------
    jplde: str
        The JPL development ephemeris version.
    """

    try:
        import lalpulsar
    except ImportError:
        raise ImportError("You must have lalsuite/lalpulsar installed to use this function")
    
    earth = body_ephemeris_path(body="earth", jplde=jplde, string=True)
    sun = body_ephemeris_path(body="sun", jplde=jplde, string=True)

    return lalpulsar.InitBarycenter(earth, sun)
