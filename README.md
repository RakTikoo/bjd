# bjd

**Batch-convert a column of observation times in a CSV or FITS table to BJD_TDB — explicitly, or not at all.**

Exoplanet and transit work lives and dies by time standards. A column of
Julian Dates is meaningless until you know its format (JD? MJD? ISO?) and
its scale (UTC? TT? TAI? TDB?) — and the wrong assumption does not crash,
it silently shifts every time by seconds to over a minute. That is enough
to fake a transit timing variation (Eastman, Siverd & Gaudi 2010, PASP,
122, 935). Existing tools make you write the Astropy incantation yourself,
every time. `bjd` is the file-level tool: point it at a table, state the
format and scale out loud, and get the same table back with a `bjd_tdb`
column appended — plus a summary showing how large the correction was, so
you can see the shift instead of trusting it.

If you do not state `--format` and `--scale`, `bjd` refuses and tells you
why. That refusal is the product.

## Demo

```console
$ bjd convert examples/times.csv --time-col jd_utc --format jd --scale utc \
      --ra 84.291188 --dec -20.273999 --site "-24.6272,-70.4042,2635" \
      --out examples/times_bjd.csv
converted 5 rows from examples/times.csv
  time column : jd_utc (jd, UTC) -> bjd_tdb (BJD_TDB)
  target      : RA 84.291188 deg, Dec -20.273999 deg (ICRS)
  site        : lat -24.627200, lon -70.404200, height 2635 m
  median shift applied: +171.94 s (min +165.99 s, max +172.06 s)
  output      : examples/times_bjd.csv
```

The output file is the input with one column appended:

```csv
jd_utc,exposure_s,note,bjd_tdb
2460000.5,60,transit ingress,2460000.5019914582
2460000.51,60,transit mid,2460000.5119907637
...
```

Forget the scale, and this happens instead:

```console
$ bjd convert examples/times.csv --time-col jd_utc --format jd \
      --ra 84.291188 --dec -20.273999 --site "-24.6272,-70.4042,2635"
error: --scale is required, and bjd refuses to guess it.
UTC, TT, TAI and TDB readings of the same instant differ by more than a
minute. Assuming the wrong scale silently shifts every time -- enough to
fake an exoplanet transit timing variation (Eastman et al. 2010). State
it explicitly: --scale utc, --scale tt, --scale tai or --scale tdb.
```

## Install

```console
$ pip install git+https://github.com/RakTikoo/bjd
```

Requires Python ≥ 3.10 and Astropy ≥ 6. Everything runs **offline**: the
conversion uses Astropy's built-in ephemeris and bundled leap-second
tables, and the tool disables network downloads for its data files.

## Usage

```
bjd convert FILE --time-col COL --format jd|mjd|iso --scale utc|tt|tai|tdb
            --ra DEG --dec DEG --site "lat,lon,height_m" [--out FILE]
bjd check
```

- `FILE` — a CSV or FITS table (Astropy reads both). All columns are
  carried through; `bjd_tdb` is appended. Output format follows the
  `--out` suffix: `.csv` for CSV, `.fits` for FITS. Without `--out`, the
  output is `<input-name>_bjd.csv` next to the input.
- `--time-col` — the column holding the observation times.
- `--format` — `jd`, `mjd`, or `iso` (ISO date-time strings such as
  `2023-02-25 00:00:00` or `2023-02-25T00:00:00`; a trailing `Z` is
  refused — state the scale with `--scale` instead).
  **Required; never guessed.**
- `--scale` — `utc`, `tt`, `tai`, or `tdb`. **Required; never guessed.**
- `--ra`, `--dec` — target coordinates in degrees, ICRS. The barycentric
  correction is a projection onto the target direction, so there is no
  sensible default.
- `--site` — the observatory as `"lat,lon,height_m"` (degrees, degrees,
  metres; height optional). Named observatories are deliberately **not**
  accepted: Astropy's site registry is a network download, and this tool
  works offline — give coordinates instead. (Paranal, used above, is
  `"-24.6272,-70.4042,2635"`. Note the order: latitude first.)

### Verify your install

```
$ bjd check
bjd check -- built-in reference case
  input   : JD 2460000.5 UTC, RA 84.291188, Dec -20.273999, site lat -24.6272, lon -70.4042, height 2635 m (Paranal)
  computed: BJD_TDB 2460000.5019914582
  expected: BJD_TDB 2460000.5019914582 (Astropy builtin ephemeris, independently reproduced by the geometric formula of Eastman et al. 2010 -- see bjd/core.py)
  result  : PASS (difference 0.0000 s, tolerance 0.1 s)
```

The expected value was computed with Astropy's `light_travel_time`
(built-in ephemeris) and independently reproduced to better than 1 ms by
the direct geometric formula of Eastman et al. (2010),
LTT = (r_obs · n̂)/c — two unrelated code paths agreeing. The derivation
is documented in `bjd/core.py`.

## Method

BJD_TDB = TDB + light-travel time, evaluated with Astropy's
`Time.light_travel_time` (built-in ephemeris): the input times are first
re-expressed in TDB, then shifted by the light-travel time from the
observatory to the solar-system barycentre along the target direction —
the definition of Eastman, Siverd & Gaudi (2010). The summary's "median
shift applied" is the difference between the output BJD_TDB values and
your input readings, so it includes the scale offset (e.g. TDB−UTC ≈
69.2 s in 2023) as well as the light-travel term (at most ≈ 500 s).

## Honest limits

- **Built-in ephemeris precision.** Astropy's built-in ephemeris is an
  analytic approximation to JPL's; Earth's barycentric position is good
  to a level that keeps the light-travel time accurate to well under
  0.1 s — far below the ~1-minute errors this tool exists to prevent.
  If you are doing sub-second transit timing for publication, cross-check
  critical epochs against a JPL-ephemeris tool such as
  [barycorrpy](https://github.com/shbhuk/barycorrpy).
- **One target, one site per run.** A table is assumed to be observations
  of a single target from a single observatory — true for essentially
  all transit follow-up tables. Mixed-target tables need one run each.
- **Mid-exposure times are your responsibility.** `bjd` converts the
  times you give it; if your column holds exposure *start* times, add
  half the exposure before converting.
- **No relativistic (Shapiro) delay or proper-motion correction.** At
  the precision of the built-in ephemeris these are below a millisecond
  for typical targets; Eastman et al. (2010) discuss when they matter.
- **Coordinates are ICRS degrees.** Precession/nutation of the target
  direction is handled inside Astropy's frame machinery; proper motion
  of the target itself is not applied.
- **Modern dates only, for the precision above.** The built-in
  ephemeris and the bundled leap-second tables are fitted/tabulated
  around the present era. Outside roughly 1900–2100 the tool still
  converts but prints a warning, and the sub-0.1 s precision claim does
  not apply (UTC itself is only defined from 1972; earlier readings are
  treated proleptically).
- **A swapped `--site` that stays in range cannot be caught.** If both
  numbers are plausible as either latitude or longitude (e.g.
  `-70.4042,-24.6272` for Paranal), the conversion succeeds with the
  wrong site and a wrong — but plausible-looking — correction. Only
  out-of-range swaps are refused. Check the printed site line.

## Tests

```console
$ python3 run_tests.py
```

26 tests (stdlib `unittest`): the pinned known-value case, both refusal
rules, JD/MJD/ISO/TT input paths, CSV and FITS round trips, the `check`
subcommand, adversarial inputs (empty file, missing column,
malformed site, site name, impossible latitude) — all of which must fail
with clean one-line errors, never tracebacks — plus the regression
tests from the 2026-11-17 adversarial review in `REVIEW.md` (non-finite
or implausible site height, ISO `T` separator, BOM-prefixed CSV,
`--out` equal to the input, one-line corrupt-file errors).

## Use with an AI agent

This repo ships a [SKILL.md](SKILL.md) written to the Agent Skills open
standard. If your coding agent supports skills, point it at this repo and
it can install bjd, run it on your tables, and read the results for you.
The skill works as documented when the repo is cloned under its real
name (`bjd`) — the standard matches a skill to its folder name.

## License

MIT. Free and open source, for science.
