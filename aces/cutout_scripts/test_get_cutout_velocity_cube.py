import sys
import types

import pytest
from astropy import constants as const
from astropy import units as u

import get_cutout_velocity_cube as line_cutout


def test_frequency_limits_use_radio_velocity_and_are_sorted():
    result = line_cutout.frequency_limits(100.0, 75.0)
    velocities = (-75.0, 75.0) * u.km / u.s
    expected = sorted((100.0 * u.GHz * (1 - velocities / const.c)).to_value(u.GHz))

    assert result == pytest.approx(expected)


@pytest.mark.parametrize("value", ["0", "-1"])
def test_positive_number_rejects_non_positive_values(value):
    with pytest.raises(Exception, match="greater than zero"):
        line_cutout.positive_number(value)


def test_wrapper_downloads_then_converts_and_overwrites(monkeypatch):
    calls = {}

    def fake_get_cutout(*args, **kwargs):
        calls["download_args"] = args
        calls["download_kwargs"] = kwargs

    class FakeCube:
        def __init__(self):
            self.allow_huge_operations = False

        def with_spectral_unit(self, unit, velocity_convention, rest_value):
            calls["large_before_velocity"] = self.allow_huge_operations
            calls["spectral_unit"] = unit
            calls["velocity_convention"] = velocity_convention
            calls["rest_value"] = rest_value
            return self

        def to(self, unit):
            calls["large_before_kelvin"] = self.allow_huge_operations
            calls["brightness_unit"] = unit
            return self

        def write(self, filename, overwrite):
            calls["write"] = (filename, overwrite)

    cube = FakeCube()
    fake_spectral_cube = types.SimpleNamespace(
        SpectralCube=types.SimpleNamespace(
            read=lambda filename: calls.update(read=filename) or cube
        )
    )

    monkeypatch.setattr(line_cutout.get_cutout, "get_cutout", fake_get_cutout)
    monkeypatch.setitem(sys.modules, "spectral_cube", fake_spectral_cube)

    line_cutout.get_cutout_velocity_cube(
        "dataset.fits", 1.2, -0.3, 20.0, 100.0, 75.0, "line.fits"
    )

    fmin, fmax = line_cutout.frequency_limits(100.0, 75.0)
    assert calls["download_args"] == ("dataset.fits", 1.2, -0.3, 20.0)
    assert calls["download_kwargs"] == {
        "fmin": pytest.approx(fmin),
        "fmax": pytest.approx(fmax),
        "outfile": "line.fits",
    }
    assert calls["read"] == "line.fits"
    assert calls["large_before_velocity"] is True
    assert calls["large_before_kelvin"] is True
    assert calls["spectral_unit"] == u.km / u.s
    assert calls["velocity_convention"] == "radio"
    assert calls["rest_value"] == 100.0 * u.GHz
    assert calls["brightness_unit"] == u.K
    assert calls["write"] == ("line.fits", True)


def test_cli_passes_all_arguments_to_wrapper(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        line_cutout,
        "get_cutout_velocity_cube",
        lambda *args: captured.update(args=args),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "get_cutout_velocity_cube.py",
            "dataset.fits",
            "1.2",
            "-0.3",
            "20",
            "100.0",
            "75",
            "line.fits",
        ],
    )

    line_cutout.main()

    assert captured["args"] == (
        "dataset.fits",
        1.2,
        -0.3,
        20.0,
        100.0,
        75.0,
        "line.fits",
    )
