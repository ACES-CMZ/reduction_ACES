import sys

import pytest

import get_cutout as cutouts


@pytest.mark.parametrize(
    ("extra", "expected"),
    [
        ([], {"vmin": None, "vmax": None, "outfile": None}),
        (["output.fits"], {"vmin": None, "vmax": None, "outfile": "output.fits"}),
        (["-10", "20"], {"vmin": -10.0, "vmax": 20.0, "outfile": None}),
        (
            ["-10", "20", "output.fits"],
            {"vmin": -10.0, "vmax": 20.0, "outfile": "output.fits"},
        ),
    ],
)
def test_existing_cli_forms_are_unchanged(monkeypatch, extra, expected):
    captured = {}
    monkeypatch.setattr(cutouts, "get_cutout", lambda *args, **kwargs: captured.update(kwargs))
    monkeypatch.setattr(
        sys,
        "argv",
        ["get_cutout.py", "dataset.fits", "1", "2", "3", *extra],
    )

    cutouts.main()

    assert captured["vmin"] == expected["vmin"]
    assert captured["vmax"] == expected["vmax"]
    assert captured["outfile"] == expected["outfile"]
    assert captured["fmin"] is None
    assert captured["fmax"] is None


@pytest.mark.parametrize("outfile", [None, "output.fits"])
def test_frequency_cli_forms(monkeypatch, outfile):
    captured = {}
    extra = ["85.5", "87.5"]
    if outfile is not None:
        extra.append(outfile)
    monkeypatch.setattr(cutouts, "get_cutout", lambda *args, **kwargs: captured.update(kwargs))
    monkeypatch.setattr(
        sys,
        "argv",
        ["get_cutout.py", "dataset.fits", "1", "2", "3", *extra, "--freq"],
    )

    cutouts.main()

    assert captured["vmin"] is None
    assert captured["vmax"] is None
    assert captured["fmin"] == 85.5
    assert captured["fmax"] == 87.5
    assert captured["outfile"] == outfile


def test_frequency_to_band_converts_ghz_and_sorts_reversed_bounds():
    from astropy import constants as const
    from astropy import units as u

    band = cutouts.frequency_to_band(100, 90)

    assert band == pytest.approx(
        sorted(
            [
                (const.c / (100 * u.GHz)).to_value(u.m),
                (const.c / (90 * u.GHz)).to_value(u.m),
            ]
        )
    )


@pytest.mark.parametrize(
    "limits", [(0, 90), (-1, 90), (float("nan"), 90), ("invalid", 90)]
)
def test_frequency_to_band_rejects_invalid_bounds(limits):
    with pytest.raises(ValueError, match="frequency limits must be"):
        cutouts.frequency_to_band(*limits)


def test_frequency_cutout_does_not_fetch_header(monkeypatch, tmp_path):
    class Response:
        url = "https://example.test/soda"
        status_code = 200
        headers = {"content-type": "application/fits"}
        content = b"FITS"

    request_params = {}

    def fake_get(url, params, timeout):
        request_params.update(params)
        return Response()

    def unexpected_header_fetch(*args, **kwargs):
        raise AssertionError("frequency mode must not fetch the FITS header")

    import requests

    monkeypatch.setattr(cutouts, "galactic_to_icrs", lambda glon, glat: (1.0, 2.0))
    monkeypatch.setattr(cutouts, "fetch_primary_header", unexpected_header_fetch)
    monkeypatch.setattr(requests, "get", fake_get)
    outfile = tmp_path / "frequency.fits"

    cutouts.get_cutout(
        "dataset.fits", 1, 2, 3, fmin=85.5, fmax=87.5, outfile=outfile
    )

    assert "BAND" in request_params
    assert outfile.read_bytes() == b"FITS"


def test_frequency_output_name_contains_ghz_range():
    name = cutouts.output_name("dataset.fits", 1, 2, 3, fmin=85.5, fmax=87.5)

    assert name == "dataset_l1_b2_r3arcsec_f85.5_87.5GHz_cutout.fits"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"fmin": 85.5},
        {"fmax": 87.5},
        {"vmin": -10, "vmax": 20, "fmin": 85.5, "fmax": 87.5},
    ],
)
def test_python_api_rejects_incomplete_or_mixed_limits(kwargs):
    with pytest.raises(ValueError):
        cutouts.get_cutout("dataset.fits", 1, 2, 3, **kwargs)
