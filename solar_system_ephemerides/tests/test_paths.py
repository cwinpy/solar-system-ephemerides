import pytest

from pathlib import Path

from solar_system_ephemerides.paths import (
    body_ephemeris_path,
    time_ephemeris_path,
    BODIES,
    JPLDE,
    TIMESPANS,
)


class TestBodyEphemerisPath:
    """
    Test the body_ephemeris_path function (which uses the BodyEphemerisPath class).
    """

    def test_no_arguments(self):
        """
        Test that body_ephemeris_path returns a TypeError with no arguments given.
        """

        with pytest.raises(TypeError):
            body_ephemeris_path()

    def test_body_setter(self):
        """
        Test setting of the body.
        """

        # TypeError for non-string inputs
        with pytest.raises(TypeError):
            body_ephemeris_path(234854)

        # ValueError for unrecognised body
        with pytest.raises(ValueError):
            body_ephemeris_path("Blah")

        # test an alias for the Earth
        b1 = body_ephemeris_path("terra", string=True)
        b2 = body_ephemeris_path("earth", string=True)

        assert b1 == b2

        # test and alias for the Sun
        b1 = body_ephemeris_path("zon", string=True)
        b2 = body_ephemeris_path("sun", string=True)

    def test_jplde_setter(self):
        """
        Test setting of JPL development ephemeris version.
        """

        with pytest.raises(ValueError):
            body_ephemeris_path("Earth", jplde=-683)

        with pytest.raises(ValueError):
            body_ephemeris_path("SUN", jplde="DE9999")

        # test passing int versus string versus DE string
        b1 = body_ephemeris_path("SUN", jplde="405", string=True)
        b2 = body_ephemeris_path("sol", jplde="DE405", string=True)
        b3 = body_ephemeris_path("sonne", jplde=405, string=True)

        assert b1 == b2 and b1 == b3

    def test_str_versus_path(self):
        """
        Check that function returns a string or a Path as expected.
        """

        strout = body_ephemeris_path("erde", jplde="200", string=True)
        pathout = body_ephemeris_path("aarde", jplde=200)

        assert isinstance(strout, str)
        assert isinstance(pathout, Path)

        assert strout == str(pathout)
        assert pathout.is_file()

    def test_timespans_setter(self):
        """
        Test setting the ephemeris timespan.
        """

        with pytest.raises(TypeError):
            body_ephemeris_path("earth", timespan=1)

        with pytest.raises(ValueError):
            body_ephemeris_path("earth", timespan="kgsdg")

        # test equivalent inputs
        b1 = body_ephemeris_path("earth", timespan="00-40", string=True)
        b2 = body_ephemeris_path("terra", timespan=(2000, 2040), string=True)
        b3 = body_ephemeris_path("ziemia", timespan=["00", 40], string=True)

        assert b1 == b2 and b1 == b3

    def test_relative_path_setter(self):
        """
        Test setting relative path.
        """

        with pytest.raises(ValueError):
            body_ephemeris_path("sol", relative_path="blah.txt")

        path = Path(__file__).parent
        b = body_ephemeris_path("sol", relative_path=path)
        assert (path / b).is_file()

    def test_all_exist(self):
        """
        Check that all the JPL ephemeris versions for all bodies for all
        timespans are actually present.
        """

        for jplde in JPLDE:
            for body in BODIES:
                for timespan in TIMESPANS:
                    b = body_ephemeris_path(body, jplde=jplde, timespan=timespan)

                    assert b.is_file()
