"""Core conversion logic for bjd.

Everything here is explicit: the caller must state the input time format
(jd / mjd / iso) and time scale (utc / tt / tai / tdb), the target
coordinates, and the observatory location. Nothing is guessed, because a
wrong guess silently shifts times by seconds to over a minute -- enough to
fake or erase an exoplanet transit timing signal.

The conversion follows Eastman, Siverd & Gaudi (2010, PASP, 122, 935),
"Achieving Better Than 1 Minute Accuracy in the Heliocentric and
Barycentric Julian Dates":  BJD_TDB = TDB + light-travel time from the
observatory to the solar-system barycentre along the target direction.
It is evaluated with Astropy's ``Time.light_travel_time`` using the
built-in ephemeris, so it works fully offline.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

# Never try to download IERS data: the tool must work offline. The bundled
# tables are more than accurate enough for this conversion (UTC -> TDB does
# not depend on UT1 at all).
iers.conf.auto_download = False

EPHEMERIS = "builtin"

FORMATS = ("jd", "mjd", "iso")
SCALES = ("utc", "tt", "tai", "tdb")

# ---------------------------------------------------------------------------
# Known-value reference case (used by `bjd check` and by the test suite).
#
# Input : JD 2460000.5 UTC  (2023-02-25 00:00:00 UTC)
# Target: RA 84.291188 deg, Dec -20.273999 deg (ICRS)
# Site  : Cerro Paranal -- lat -24.6272, lon -70.4042, height 2635 m
#
# Expected BJD_TDB = 2460000.5019914582
#   (light-travel time +102.876705 s; total shift vs the JD_UTC reading
#    +172.062 s, of which 69.185 s is TDB-UTC.)
#
# How the expected value was derived: computed with Astropy 8.0.2
# (Time.light_travel_time, ephemeris='builtin') AND independently
# reproduced to better than 1 ms by the direct geometric formula of
# Eastman et al. (2010):  LTT = (r_obs . n_hat) / c, with r_obs the
# observatory's barycentric position and n_hat the unit vector to the
# target. Two unrelated code paths agreeing to <1 ms is the reference.
# ---------------------------------------------------------------------------
KNOWN_INPUT_JD = 2460000.5
KNOWN_INPUT_SCALE = "utc"
KNOWN_RA = 84.291188
KNOWN_DEC = -20.273999
KNOWN_SITE = (-24.6272, -70.4042, 2635.0)  # lat deg, lon deg, height m
KNOWN_BJD_TDB = 2460000.5019914582
KNOWN_TOLERANCE_S = 0.1


class BjdError(Exception):
    """A user-facing error: print the message, exit non-zero, no traceback."""


def parse_site(text):
    """Parse a --site value into an EarthLocation.

    Accepts "lat,lon" or "lat,lon,height_m" (degrees, degrees, metres).
    Named observatories are NOT accepted: Astropy's site registry needs a
    network download, and this tool works offline -- it says so instead of
    failing obscurely.
    """
    if text is None:
        raise BjdError("--site is required: give the observatory as "
                       '"lat,lon,height_m" (degrees, degrees, metres).')
    parts = [p.strip() for p in str(text).split(",")]
    if len(parts) not in (2, 3) or any(p == "" for p in parts):
        # A single token is probably a site name -- explain, don't crash.
        if len(parts) == 1 and not _looks_like_number(parts[0]):
            raise BjdError(
                f'--site "{text}" looks like a name. Named observatories '
                "need a network download in Astropy, and bjd works offline. "
                'Please give coordinates instead: --site "lat,lon,height_m" '
                '(for example Paranal is "-24.6272,-70.4042,2635").')
        raise BjdError(f'--site "{text}" is malformed: expected '
                       '"lat,lon,height_m" (for example '
                       '"-24.6272,-70.4042,2635").')
    try:
        lat = float(parts[0])
        lon = float(parts[1])
        height = float(parts[2]) if len(parts) == 3 else 0.0
    except ValueError:
        raise BjdError(f'--site "{text}" is malformed: lat, lon and height '
                       "must be numbers.")
    if not -90.0 <= lat <= 90.0:
        raise BjdError(f"--site latitude {lat} is outside [-90, 90] degrees. "
                       'Remember the order is "lat,lon,height_m".')
    if not -180.0 <= lon <= 180.0:
        raise BjdError(f"--site longitude {lon} is outside [-180, 180] "
                       "degrees.")
    # Height must be finite and plausible. A nan/inf height used to pass
    # through EarthLocation and silently produce an all-nan bjd_tdb
    # column with exit code 0 (found in adversarial review, 2026-11-17);
    # a height of e.g. 1e9 m silently shifted results by seconds. The
    # range below covers every ground observatory (Dead Sea shore
    # ~-430 m to Chajnantor ~5600 m), aircraft and stratospheric
    # balloons, while refusing values that cannot be a real site.
    if not np.isfinite(height):
        raise BjdError(f"--site height {parts[2] if len(parts) == 3 else height} "
                       "is not a finite number of metres.")
    if not -10000.0 <= height <= 100000.0:
        raise BjdError(f"--site height {height:g} m is implausible for an "
                       "observatory (allowed range -10000 to 100000 m). "
                       'Remember the order is "lat,lon,height_m", with '
                       "height in metres.")
    return EarthLocation(lat=lat * u.deg, lon=lon * u.deg,
                         height=height * u.m)


def _looks_like_number(text):
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def validate_target(ra, dec):
    if ra is None or dec is None:
        raise BjdError("--ra and --dec are required: the barycentric "
                       "correction depends on the target direction.")
    try:
        ra = float(ra)
        dec = float(dec)
    except (TypeError, ValueError):
        raise BjdError("--ra and --dec must be numbers (degrees, ICRS).")
    if not np.isfinite(ra) or not 0.0 <= ra < 360.0:
        raise BjdError(f"--ra {ra} is outside [0, 360) degrees.")
    if not np.isfinite(dec) or not -90.0 <= dec <= 90.0:
        raise BjdError(f"--dec {dec} is outside [-90, 90] degrees.")


def to_bjd_tdb(values, fmt, scale, ra_deg, dec_deg, location):
    """Convert an array of times to BJD_TDB.

    ``values`` is a sequence of floats for fmt jd/mjd, or strings for iso.
    Returns (bjd_times, input_times) as Astropy Time objects.
    """
    # Validate everything here as well as in the CLI, so direct library
    # callers get the same clean BjdError refusals instead of raw Astropy
    # exceptions -- or worse, a silent nan result for a nan RA.
    if fmt not in FORMATS:
        raise BjdError(f"format {fmt!r} is not one of {FORMATS}: state the "
                       "time format explicitly (jd, mjd or iso).")
    if scale not in SCALES:
        raise BjdError(f"scale {scale!r} is not one of {SCALES}: state the "
                       "time scale explicitly (utc, tt, tai or tdb).")
    validate_target(ra_deg, dec_deg)
    # Masked (missing) elements must be refused before float() touches
    # them: float(masked) leaks a UserWarning to stderr on the way to
    # the (correct) non-finite error below.
    mask = getattr(values, "mask", None)
    if mask is not None:
        try:
            has_masked = bool(np.any(mask))
        except Exception:
            has_masked = False
        if has_masked:
            raise BjdError("the time column contains missing values; "
                           "fix or remove those rows first.")
    if fmt in ("jd", "mjd"):
        try:
            numeric = np.asarray([float(v) for v in values], dtype=float)
        except (TypeError, ValueError):
            raise BjdError("the time column contains values that are not "
                           f"numbers, but --format {fmt} needs numbers.")
        if numeric.size == 0:
            raise BjdError("the time column is empty: nothing to convert.")
        if not np.all(np.isfinite(numeric)):
            raise BjdError("the time column contains missing or non-finite "
                           "values; fix or remove those rows first.")
        times = Time(numeric, format=fmt, scale=scale)
    else:  # iso
        strings = [str(v).strip() for v in values]
        if not strings or any(s == "" for s in strings):
            raise BjdError("the time column contains empty values; "
                           "--format iso needs a date-time in every row.")
        # Accept both ISO 8601 separators: "2023-02-25 00:00:00" (Astropy
        # 'iso') and "2023-02-25T00:00:00" (Astropy 'isot') -- the T form
        # is what most pipelines emit, and refusing it served no purpose.
        # A trailing 'Z' is refused explicitly (rather than depending on
        # Astropy's isot parser, which accepts it only for UTC): the
        # scale is stated with --scale, and a redundant zone designator
        # adds a second place for the two to disagree.
        if any(s.endswith(("Z", "z")) for s in strings):
            raise BjdError(
                "the time column contains ISO date-times ending in 'Z'; "
                "remove the 'Z' and state the scale with --scale "
                "('--scale utc' if the times are UTC).")
        try:
            times = Time(strings, format="iso", scale=scale)
        except Exception:
            try:
                times = Time(strings, format="isot", scale=scale)
            except Exception:
                raise BjdError(
                    "could not parse the time column as ISO date-times; "
                    "expected forms like '2023-02-25 00:00:00' or "
                    "'2023-02-25T00:00:00' (no trailing 'Z' -- state the "
                    "scale with --scale).")
    target = SkyCoord(ra=float(ra_deg) * u.deg, dec=float(dec_deg) * u.deg,
                      frame="icrs")
    ltt = times.light_travel_time(target, location=location,
                                  ephemeris=EPHEMERIS)
    bjd = times.tdb + ltt
    if not np.all(np.isfinite(bjd.jd)):
        raise BjdError("the conversion produced non-finite BJD_TDB "
                       "values; check the site coordinates and time "
                       "values for missing or invalid entries.")
    return bjd, times


def run_known_case():
    """Convert the built-in reference case; return (computed, expected)."""
    location = EarthLocation(lat=KNOWN_SITE[0] * u.deg,
                             lon=KNOWN_SITE[1] * u.deg,
                             height=KNOWN_SITE[2] * u.m)
    bjd, _ = to_bjd_tdb([KNOWN_INPUT_JD], "jd", KNOWN_INPUT_SCALE,
                        KNOWN_RA, KNOWN_DEC, location)
    return float(bjd.jd[0]), KNOWN_BJD_TDB
