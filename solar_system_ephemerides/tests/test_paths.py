import pytest

from pathlib import Path

from solar_system_ephemerides.paths import body_ephemeris_path, time_ephemeris_path


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
