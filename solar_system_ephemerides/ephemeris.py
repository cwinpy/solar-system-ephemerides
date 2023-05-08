import numpy as np
from gzip import open as gzopen

from astropy import constants
from astropy import units as u

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

        self.pos = data[:, 1:4]  # position in lt/s
        self.vel = data[:, 4:7]  # velocity n lt/s
        self.acc = data[:, 7:10]  # acceleration in lt/s^2

    @property
    def pos(self):
        # position in lt/s
        return self.get_pos()

    @pos.setter
    def pos(self, pos):
        self._pos = np.asarray(pos, dtype=float)

        if self._pos.shape[-1] != 3:
            raise ValueError("Position must contain x, y and z components")

    def get_pos(
        self, coords: str = "xyz", si: bool = False, astropyunits: bool = False
    ):
        """
        Return a numpy array with the position at each time step. By
        default, these will be in light seconds.

        Parameters
        ----------
        coords: str
            The coordinate(s) to return. Default is "xyz", i.e., the 3D
            position.
        si: bool
            Return the position in metres rather than light seconds. Default is
            False.
        astropyunits: bool
            Return the array as an astropy Quantity with the associated units.
            Default is False.
        """

        try:
            self._pos
        except AttributeError:
            return None

        xyz = {"x": 0, "y": 1, "z": 2}
        if not (set("xyz") & set(coords)):
            raise ValueError("Coordinates must contain 'x', 'y', and/or 'z'")

        sl = np.r_[[xyz[v] for v in set(coords) if v in xyz]]
        if not si:
            return self._pos[:, sl] * u.lsec if astropyunits else self._pos[:, sl]
        else:
            return (
                self._pos[:, sl] * u.s * constants.c
                if astropyunits
                else self._pos[:, sl] * constants.c.value
            )

    @property
    def pos_si(self):
        # return position in metres
        return self.get_pos(si=True) if self.pos is not None else None

    @property
    def pos_x(self):
        # return x position
        return self.get_pos(coords="x")

    @property
    def pos_y(self):
        # return y position
        return self.get_pos(coords="y")

    @property
    def pos_z(self):
        # return z position
        return self.get_pos(coords="z")

    @property
    def vel(self):
        return self.get_vel()

    @vel.setter
    def vel(self, vel):
        self._vel = np.asarray(vel, dtype=float)

        if self._vel.shape[-1] != 3:
            raise ValueError("Velocity must contain x, y and z components")

    def get_vel(
        self, coords: str = "xyz", si: bool = False, astropyunits: bool = False
    ):
        """
        Return a numpy array with the velocity at each time step. By
        default, these will be in light seconds/second.

        Parameters
        ----------
        coords: str
            The coordinate(s) to return. Default is "xyz", i.e., the 3D
            velocity.
        si: bool
            Return the velocity in metres/s rather than light seconds/s.
            Default is False.
        astropyunits: bool
            Return the array as an astropy Quantity with the associated units.
            Default is False.
        """

        try:
            self._vel
        except AttributeError:
            return None

        xyz = {"x": 0, "y": 1, "z": 2}
        if not (set("xyz") & set(coords)):
            raise ValueError("Coordinates must contain 'x', 'y', and/or 'z'")

        sl = np.r_[[xyz[v] for v in set(coords) if v in xyz]]
        if not si:
            return self._vel[:, sl] * u.lsec / u.s if astropyunits else self._vel[:, sl]
        else:
            return (
                self._vel[:, sl] * constants.c
                if astropyunits
                else self._vel[:, sl] * constants.c.value
            )

    @property
    def vel_si(self):
        # return velocity in metres/s
        return self.get_vel(si=True)

    @property
    def vel_x(self):
        # return x velocity
        return self.get_vel(coords="x")

    @property
    def vel_y(self):
        # return y velocity
        return self.get_vel(coords="y")

    @property
    def vel_z(self):
        # return z velocity
        return self.get_vel(coords="z")

    @property
    def acc(self):
        return self.get_acc()

    @acc.setter
    def acc(self, acc):
        self._acc = np.asarray(acc, dtype=float)

        if self._acc.shape[-1] != 3:
            raise ValueError("Acceleration must contain x, y and z components")

    def get_acc(
        self, coords: str = "xyz", si: bool = False, astropyunits: bool = False
    ):
        """
        Return a numpy array with the acceleration at each time step. By
        default, these will be in light seconds/second.

        Parameters
        ----------
        coords: str
            The coordinate(s) to return. Default is "xyz", i.e., the 3D
            acceleration.
        si: bool
            Return the acceleration in metres/s^2 rather than light
            seconds/s^2. Default is False.
        astropyunits: bool
            Return the array as an astropy Quantity with the associated units.
            Default is False.
        """

        try:
            self._acc
        except AttributeError:
            return None

        xyz = {"x": 0, "y": 1, "z": 2}
        if not (set("xyz") & set(coords)):
            raise ValueError("Coordinates must contain 'x', 'y', and/or 'z'")

        sl = np.r_[[xyz[v] for v in set(coords) if v in xyz]]
        if not si:
            return (
                self._acc[:, sl] * u.lsec / u.s**2
                if astropyunits
                else self._acc[:, sl]
            )
        else:
            return (
                self._acc[:, sl] * constants.c / u.s
                if astropyunits
                else self._acc[:, sl] * constants.c.value
            )

    @property
    def acc_si(self):
        # return acceleration in metres/s^2
        return self.get_acc(si=True)

    @property
    def acc_x(self):
        # return x acceleration
        return self.get_acc(coords="x")

    @property
    def acc_y(self):
        # return y acceleration
        return self.get_acc(coords="y")

    @property
    def acc_z(self):
        # return z acceleration
        return self.get_acc(coords="z")

    @property
    def times(self):
        try:
            return np.arange(
                self.t0, self.t0 + self.timestep * self.nentries, self.nentries
            )
        except AttributeError:
            return None

    def write(self, outfile, header=None):
        """
        Write the ephemeris out to a file. If the file extension is ".txt" or
        ".dat" this will be output to a plain ascii file. If it ends with ".gz"
        then file will be gzipped.

        Parameters
        ----------
        outfile: str
            The path to the output file.
        header: str
            A string to output at the header in the file. Each line in this
            should start with a "#" to denote a comment line.
        """

        if outfile.endswith(".hdf5"):
            raise NotImplementedError("Output to HDF5 files is not yet implemented.")

        # gzip if extension ends with '.gz'
        if outfile.endswith(".gz"):
            try:
                import gzip

                fp = gzip.open(outfile, "wb")
            except IOError:
                Exception("Problem opening gzip output file '{}'".format(outfile))
        else:
            try:
                fp = open(outfile, "w")
            except IOError:
                Exception("Problem opening output file '{}'".format(outfile))

        try:
            _ = self.t0
            _ = self.pos
        except AttributeError:
            raise RuntimeError("Cannot write ephemeris as it has not been filled in yet.")

        # write out header
        fp.write(header + ("\n" if not header.endswith("\n") else ""))

        # write out start time (GPS), time interval (secs), number of entries
        fp.write(f"{int(self.t0)}\t{self.timestep}\t{len(self)}\n")

        times = self.times
        pos = self.pos
        vel = self.vel
        acc = self.acc
        for i in range(len(self)):
            fp.write(f"{times[i]:.6f}\t")
            fp.write(
                f"{pos[i, 0]:.16e}\t{pos[i, 1]:.16e}\t{pos[i, 2]:.16e}\t"
            )
            fp.write(
                f"{vel[i, 0]:.16e}\t{vel[i, 1]:.16e}\t{vel[i, 2]:.16e}\t"
            )
            fp.write(
                f"{acc[i, 0]:.16e}\t{acc[i, 1]:.16e}\t{acc[i, 2]:.16e}\n"
            )

        fp.close()

    def __len__(self):
        try:
            return self.nentries
        except AttributeError:
            return 0


def lal_ephemeris_data(jplde: str = "DE405", timespan="00-40"):
    """
    Return a LAL EphemerisData object containing the ephemerides for the Earth
    and Sun.

    Parameters
    ----------
    jplde: str
        The JPL development ephemeris version. Default is DE405.
    timespan: str
        The timespan of the ephemeris. This defaults to "00-40", which is
        currently the only time span available within this package.
    """

    try:
        import lalpulsar
    except ImportError:
        raise ImportError(
            "You must have lalsuite/lalpulsar installed to use this function"
        )

    earth = body_ephemeris_path(
        body="earth", jplde=jplde, timespan=timespan, string=True
    )
    sun = body_ephemeris_path(body="sun", jplde=jplde, timespan=timespan, string=True)

    return lalpulsar.InitBarycenter(earth, sun)


def lal_time_ephemeris_data(units: str = "TCB"):
    """
    Return a LAL TimeCorrectionData object contain the time correction data.

    Parameters
    ----------
    units: str
        The time system, either TCB (the default) or TDB.
    """

    try:
        import lalpulsar
    except ImportError:
        raise ImportError(
            "You must have lalsuite/lalpulsar installed to use this function"
        )

    time = time_ephemeris_path(units=units, string=True)

    return lalpulsar.InitTimeCorrections(time)
