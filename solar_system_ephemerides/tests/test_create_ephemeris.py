"""
Test the creation of an ephemeris.

The DE200, DE405, DE421 and DE430 ephemerides in this repository are those that
were originally contained within LALSuite. We can therefore test that this code
reproduces a couple of them.  
"""

from pathlib import Path
import shutil

import pytest
import numpy as np

from solar_system_ephemerides.ephemeris import (
    BodyEphemeris,
    lal_ephemeris_data,
    lal_time_ephemeris_data,
)
from solar_system_ephemerides.generate import generate_ephemeris
from solar_system_ephemerides.paths import time_ephemeris_path


class TestGenerateEphemeris:
    @classmethod
    def setup_class(cls):
        # create output directory
        cls.outputdir = Path("test_ephemerides")
        cls.outputdir.mkdir(parents=True, exist_ok=True)

    @classmethod
    def teardown_class(cls):
        """
        Remove test output plots.
        """

        shutil.rmtree(cls.outputdir)

    def test_exceptions(self):
        """
        Test various exceptions within the function.
        """

        with pytest.raises(ValueError):
            generate_ephemeris()

        with pytest.raises(ValueError):
            generate_ephemeris(jplde="DE999")

        with pytest.raises(ValueError):
            generate_ephemeris(jplde="DE405", body="lsgdglfad")

        with pytest.raises(ValueError):
            generate_ephemeris(
                jplde="DE405", body="sun", yearstart=2000, gpsstart=900000000
            )

        with pytest.raises(ValueError):
            generate_ephemeris(jplde="DE405", body="sun", gpsstart="kfsfd")

        with pytest.raises(ValueError):
            generate_ephemeris(jplde="DE405", body="sun", yearstart="kfsfd")

        with pytest.raises(ValueError):
            generate_ephemeris(
                jplde="DE405", body="sun", yearstart=2000, interval="ksgdkc"
            )

        with pytest.raises(ValueError):
            generate_ephemeris(
                jplde="DE405", body="sun", yearstart=2000, interval=8, nyears="kslfdfvc"
            )

    def test_generate_de200_earth(self):
        """
        Generate a DE200 ephemeris for the Earth and compare with the one in
        the package.
        """

        outfile = self.outputdir / "earth00-40-DE200.dat.gz"

        ephem = generate_ephemeris(
            body="earth",
            jplde="DE200",
            nyears=40,
            interval=2,
            yearstart=2000,
            output=outfile,
        )

        assert isinstance(ephem, BodyEphemeris)

        # read in the packaged file
        ephempackaged = BodyEphemeris(body="earth", jplde="DE200", timespan="00-40")

        assert np.allclose(ephem.times, ephempackaged.times, atol=1e-12)
        assert np.allclose(ephem.pos, ephempackaged.pos, atol=1e-12)
        assert np.allclose(ephem.vel, ephempackaged.vel, atol=1e-12)
        assert np.allclose(ephem.acc, ephempackaged.acc, atol=1e-12)

        # read in output file and recheck
        ephemreread = BodyEphemeris(path=outfile)

        assert np.allclose(ephemreread.times, ephempackaged.times, atol=1e-12)
        assert np.allclose(ephemreread.pos, ephempackaged.pos, atol=1e-12)
        assert np.allclose(ephemreread.vel, ephempackaged.vel, atol=1e-12)
        assert np.allclose(ephemreread.acc, ephempackaged.acc, atol=1e-12)

    def test_generate_de405_sun(self):
        """
        Generate a DE405 ephemeris for the Sun and compare with the one in
        the package.
        """

        outfile = self.outputdir / "sun00-40-DE405.dat.gz"

        ephem = generate_ephemeris(
            body="sun",
            jplde="DE405",
            nyears=40,
            interval=20,
            yearstart=2000,
            output=outfile,
        )

        assert isinstance(ephem, BodyEphemeris)

        # read in the packaged file
        ephempackaged = BodyEphemeris(body="sun", jplde="DE405", timespan="00-40")

        assert np.allclose(ephem.times, ephempackaged.times, atol=1e-12)
        assert np.allclose(ephem.pos, ephempackaged.pos, atol=1e-12)
        assert np.allclose(ephem.vel, ephempackaged.vel, atol=1e-12)
        assert np.allclose(ephem.acc, ephempackaged.acc, atol=1e-12)

        # read in output file and recheck
        ephemreread = BodyEphemeris(path=outfile)

        assert np.allclose(ephemreread.times, ephempackaged.times, atol=1e-12)
        assert np.allclose(ephemreread.pos, ephempackaged.pos, atol=1e-12)
        assert np.allclose(ephemreread.vel, ephempackaged.vel, atol=1e-12)
        assert np.allclose(ephemreread.acc, ephempackaged.acc, atol=1e-12)


def test_lal_ephemeris_data():
    """
    Test using lalpulsar to read in the data in the package.
    """

    edat = lal_ephemeris_data(jplde="DE421", timespan="00-40")

    # Earth values
    posE = np.array([e.pos for e in edat.ephemE])
    velE = np.array([e.vel for e in edat.ephemE])
    accE = np.array([e.acc for e in edat.ephemE])

    # Sun values
    posS = np.array([e.pos for e in edat.ephemS])
    velS = np.array([e.vel for e in edat.ephemS])
    accS = np.array([e.acc for e in edat.ephemS])

    earth = BodyEphemeris(body="earth", jplde="DE421", timespan="00-40")

    assert np.allclose(posE, earth.pos, atol=1e-12)
    assert np.allclose(velE, earth.vel, atol=1e-12)
    assert np.allclose(accE, earth.acc, atol=1e-12)

    sun = BodyEphemeris(body="Sun", jplde=421, timespan="00-40")

    assert np.allclose(posS, sun.pos, atol=1e-12)
    assert np.allclose(velS, sun.vel, atol=1e-12)
    assert np.allclose(accS, sun.acc, atol=1e-12)


def test_lal_time_ephemeris_data():
    tdb = lal_time_ephemeris_data("TDB")
    tcb = lal_time_ephemeris_data("TCB")

    assert tdb.timeEphemeris == time_ephemeris_path(units="TDB", string=True)
    assert tcb.timeEphemeris == time_ephemeris_path(units="TCB", string=True)
