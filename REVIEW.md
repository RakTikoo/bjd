# Adversarial review — bjd (2026-11-17)

Independent review (reviewer did not build the tool). Gauntlet:
adversarial probing, an independent oracle check, a fresh-venv install
with the README demo reproduced, and the full test battery. The gauntlet was run locally;
at review time nothing had been published or pushed. `--site` is lat,lon,height by
documented design — not changed.

## Headline verification (oracle)

BJD_TDB was recomputed from scratch with the geometric formula of
Eastman et al. (2010), LTT = (r_obs · n̂)/c, using
`get_body_barycentric_posvel` + `EarthLocation.get_gcrs_posvel` —
never `Time.light_travel_time`, the tool's own code path.

| Case | Oracle BJD_TDB | Tool BJD_TDB | Diff |
|---|---|---|---|
| JD 2460000.5 UTC, Paranal, RA 84.291188 / Dec −20.273999 | 2460000.5019914582372 | 2460000.5019914582372 | 0.00 µs |
| JD 2460600.25 UTC, Mauna Kea, RA 213.704 / Dec 12.345 | 2460600.2456050626934 | 2460600.2456050626934 | 0.00 µs |
| JD 2459000.75 UTC, SAAO, RA 30 / Dec −60 | 2459000.7505483595 | 2459000.7505483595 | <0.001 ms |
| JD 2460200.0 UTC, Greenwich, RA 0 / Dec 0 | 2460200.0064613814 | 2460200.0064613814 | <0.001 ms |

Oracle LTT for the reference case: +102.876705 s and TDB−UTC 69.1853 s,
matching the constants documented in `bjd/core.py`. The pinned known
value 2460000.5019914582 is correct. (Both paths share Astropy's
built-in ephemeris for Earth's position, so this verifies the
conversion mathematics, not the ephemeris against JPL — the README
already says to cross-check publication-grade epochs with a
JPL-ephemeris tool.)

## Test battery

- Before fixes: `run_tests.py` — **15 tests, all passed**.
- After fixes: **26 tests, all passed** (15 original + 11 regression
  tests added by this review).

## Fresh install + demo

- `pip install` into a fresh venv succeeded (astropy 8.0.2, bjd
  0.1.0). Note: the venv could not be created in /tmp — /tmp is a
  512 MB tmpfs shared with other work and the install ran out of
  space twice (`Errno 28`); a fresh venv on disk was used instead and
  removed afterwards.
- README demo reproduced **exactly**: stdout identical to the README
  block (median shift +171.94 s, min +165.99 s, max +172.06 s), output
  CSV byte-identical to the committed `examples/times_bjd.csv`, the
  missing-`--scale` refusal matches, and `bjd check` PASSes. Verified
  both before and after the fixes below.

## Findings and fixes

### HIGH — non-finite site height silently produced NaN results
`--site "-24.6,-70.4,nan"` (also `inf`/`-inf`) exited 0 and wrote a
`bjd_tdb` column of `nan`, with a garbled site line (lat 90, lon 0) —
silent corruption of exactly the kind this tool exists to prevent.
**Fixed:** `parse_site` refuses non-finite heights.

### MEDIUM — implausible site height silently shifted results
Height `1e9` m was accepted and shifted BJD_TDB by ~+3.3 s; `-100000` m
was accepted too. **Fixed:** heights outside −10000…100000 m are
refused with an explanatory error (covers every ground observatory,
aircraft and stratospheric balloons).

### MEDIUM — ISO `T` separator rejected
`2023-02-25T00:00:00` — the commonest ISO 8601 form — was refused.
**Fixed:** both `iso` and `isot` forms are accepted. A trailing `Z`
is refused explicitly with a clean one-line error (state the scale
with `--scale`); this is now documented in the README.

### MEDIUM — library API bypassed the CLI's validation
Calling `core.to_bjd_tdb` directly: an unknown format fell through to
the ISO branch with a misleading error; an unknown scale raised a raw
`ScaleValueError`; Dec 120° raised a raw `ValueError`; NaN RA silently
returned a NaN BJD. **Fixed:** `to_bjd_tdb` validates format, scale
and target itself (raising `BjdError`), `validate_target` handles
non-numeric input, and a non-finite conversion result is refused.

### LOW — corrupt/unreadable files produced multi-line error dumps
An unidentifiable file embedded Astropy's whole format registry in
the error, contradicting the README's "clean one-line errors" claim;
ISO parse errors were multi-line too. **Fixed:** read errors are
truncated to the first line and ISO errors are a single composed line.

### LOW — BOM-prefixed CSVs (Excel exports) failed
A UTF-8 BOM became part of the first column name, so a correct
`--time-col jd_utc` was reported "not found". **Fixed:** a leading BOM
is stripped from the first column name after reading.

### LOW — masked time values leaked a UserWarning to stderr
A missing time cell in an otherwise valid CSV printed a Python
`UserWarning` (with a source path/line) before the clean error.
**Fixed:** masked values are detected and refused before conversion.

### LOW — `--out` equal to the input overwrote the original table
**Fixed:** refused with a clean error; the input file is left intact.

### LOW — extreme dates converted silently
JD 0, JD −1000 and JD 9999999 were converted with only Astropy's
cryptic IERS warning; the README's precision claims do not hold there
and no date-range limit was documented. **Fixed (conservatively):**
conversion still proceeds (the dates were stated explicitly), but a
plain-language stderr warning is printed for times outside
1900–2100, and the README's Honest limits now state the modern-era
precision scope. No date is refused.

## README accuracy notes

- Demo, refusal text, `bjd check` output, method description and the
  0.1 s tolerance framing were all reproduced and are accurate.
- "clean one-line errors" was overstated for corrupt files / ISO
  errors — now true after the fixes above.
- Honest limits were missing: the modern-date precision scope and the
  fact that an in-range swapped lat/lon cannot be detected. Both
  bullets added. The `--site` order warning ("latitude first") was
  already present.
- The install line (`pip install git+https://github.com/RakTikoo/bjd`)
  could not be verified from local work alone; the identical local
  tree installs cleanly.

## Not changed (by design)

- `--site` remains lat,lon,height — documented design decision.
- Reading JD numbers as MJD (or vice versa) when the user states the
  wrong `--format` still "succeeds": that is the tool's explicit-format
  contract, not a bug.
- Date-only ISO strings (`2023-02-25`) parse as midnight, standard
  Astropy ISO behaviour.
