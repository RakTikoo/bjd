"""Command-line interface for bjd.

Subcommands:
  bjd convert FILE --time-col COL --format jd|mjd|iso --scale utc|tt|tai|tdb
              --ra DEG --dec DEG --site "lat,lon,height_m" [--out FILE]
  bjd check   -- verify the converter against a built-in known value

--format and --scale have NO defaults on purpose. If either is missing the
tool refuses, explains why, and exits: guessing a time standard can shift
results by over a minute, silently.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

from . import core
from .core import BjdError

REFUSAL_FORMAT = (
    "error: --format is required, and bjd refuses to guess it.\n"
    "The same column of numbers means different instants read as JD, MJD "
    "or ISO text (JD and MJD differ by 2,400,000.5 days). State it "
    "explicitly: --format jd, --format mjd or --format iso.")

REFUSAL_SCALE = (
    "error: --scale is required, and bjd refuses to guess it.\n"
    "UTC, TT, TAI and TDB readings of the same instant differ by more "
    "than a minute. Assuming the wrong scale silently shifts every time "
    "-- enough to fake an exoplanet transit timing variation "
    "(Eastman et al. 2010). State it explicitly: --scale utc, --scale tt, "
    "--scale tai or --scale tdb.")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="bjd",
        description="Batch-convert observation times in a CSV or FITS "
                    "table to BJD_TDB -- explicitly, or not at all.")
    sub = parser.add_subparsers(dest="command")

    conv = sub.add_parser(
        "convert", help="convert a time column in a CSV/FITS table")
    conv.add_argument("file", help="input table (CSV or FITS)")
    conv.add_argument("--time-col", required=True,
                      help="name of the column holding the observation "
                           "times")
    conv.add_argument("--format", default=None, choices=core.FORMATS,
                      help="format of the time column: jd, mjd or iso "
                           "(required -- never guessed)")
    conv.add_argument("--scale", default=None, choices=core.SCALES,
                      help="time scale of the time column: utc, tt, tai "
                           "or tdb (required -- never guessed)")
    conv.add_argument("--ra", type=float, required=True,
                      help="target right ascension, degrees (ICRS)")
    conv.add_argument("--dec", type=float, required=True,
                      help="target declination, degrees (ICRS)")
    conv.add_argument("--site", required=True,
                        help='observatory as "lat,lon,height_m" '
                             '(e.g. "-24.6272,-70.4042,2635" for Paranal)')
    conv.add_argument("--out", default=None,
                      help="output file (default: <input>_bjd.csv; use a "
                           ".fits suffix for FITS output)")

    sub.add_parser("check", help="verify the converter against a built-in "
                                 "known BJD_TDB value")
    return parser


def _read_table(path):
    p = Path(path)
    if not p.is_file():
        raise BjdError(f"input file not found: {path}")
    if p.stat().st_size == 0:
        raise BjdError(f"input file is empty: {path}")
    try:
        table = Table.read(p)
    except Exception as exc:
        # Astropy's unidentifiable-format error embeds its whole format
        # registry (dozens of lines); keep the user-facing error to the
        # first line, as the README promises clean one-line errors.
        first_line = str(exc).strip().splitlines()[0] if str(exc).strip() \
            else type(exc).__name__
        raise BjdError(f"could not read {path} as a CSV or FITS table: "
                       f"{first_line}")
    # Excel-exported CSVs often start with a UTF-8 byte-order mark, which
    # Astropy keeps as part of the first column name ("\ufeffjd_utc"),
    # so the user's correct --time-col is reported as not found. Strip it.
    if table.colnames and table.colnames[0].startswith("\ufeff"):
        table.rename_column(table.colnames[0],
                            table.colnames[0].lstrip("\ufeff"))
    if len(table) == 0:
        raise BjdError(f"{path} contains no data rows: nothing to convert.")
    return table


def cmd_convert(args):
    # The refusal rule: format and scale must be stated, never defaulted.
    if args.format is None:
        print(REFUSAL_FORMAT, file=sys.stderr)
        return 2
    if args.scale is None:
        print(REFUSAL_SCALE, file=sys.stderr)
        return 2

    core.validate_target(args.ra, args.dec)
    location = core.parse_site(args.site)
    table = _read_table(args.file)

    if args.time_col not in table.colnames:
        raise BjdError(f"column '{args.time_col}' not found in "
                       f"{args.file}. Available columns: "
                       + ", ".join(table.colnames))
    if "bjd_tdb" in table.colnames:
        raise BjdError("the input already has a 'bjd_tdb' column; "
                       "refusing to overwrite it. Rename or remove it "
                       "first.")

    bjd, times = core.to_bjd_tdb(table[args.time_col], args.format,
                                 args.scale, args.ra, args.dec, location)
    table["bjd_tdb"] = bjd.jd

    if args.out:
        out_path = Path(args.out)
    else:
        out_path = Path(args.file).with_name(
            Path(args.file).stem + "_bjd.csv")
    # Never overwrite the input file itself: the original table is the
    # user's only copy of the unconverted times in many pipelines.
    try:
        same_file = out_path.resolve() == Path(args.file).resolve()
    except OSError:
        same_file = False
    if same_file:
        raise BjdError("--out is the input file itself; refusing to "
                       "overwrite the original table. Choose a different "
                       "--out path.")
    try:
        table.write(out_path, overwrite=True)
    except Exception as exc:
        raise BjdError(f"could not write output file {out_path}: {exc}")

    # The precision claims in the README hold for the modern era: UTC
    # leap seconds are only tabulated up to the bundled table's horizon
    # and the built-in ephemeris is an analytic theory fitted around the
    # present epoch. Far outside it, still convert -- the user stated
    # the dates explicitly -- but say plainly that precision is not
    # guaranteed instead of letting Astropy's IERS warning be the only
    # (cryptic) signal.
    approx_year = 2000.0 + (times.tdb.jd - 2451545.0) / 365.25
    if np.any(approx_year < 1900.0) or np.any(approx_year > 2100.0):
        print("warning: some input times fall outside 1900-2100; the "
              "built-in ephemeris and bundled leap-second tables do not "
              "guarantee the usual precision there.", file=sys.stderr)

    shift_s = (bjd.jd - times.jd) * 86400.0
    lat = float(location.lat.deg)
    lon = float(location.lon.deg)
    height = float(location.height.to_value("m"))
    print(f"converted {len(table)} rows from {args.file}")
    print(f"  time column : {args.time_col} "
          f"({args.format}, {args.scale.upper()}) -> bjd_tdb (BJD_TDB)")
    print(f"  target      : RA {args.ra} deg, Dec {args.dec} deg (ICRS)")
    print(f"  site        : lat {lat:.6f}, lon {lon:.6f}, "
          f"height {height:.0f} m")
    print(f"  median shift applied: {np.median(shift_s):+.2f} s "
          f"(min {np.min(shift_s):+.2f} s, max {np.max(shift_s):+.2f} s)")
    print(f"  output      : {out_path}")
    return 0


def cmd_check(_args):
    computed, expected = core.run_known_case()
    diff_s = abs(computed - expected) * 86400.0
    print("bjd check -- built-in reference case")
    print(f"  input   : JD {core.KNOWN_INPUT_JD} "
          f"{core.KNOWN_INPUT_SCALE.upper()}, "
          f"RA {core.KNOWN_RA}, Dec {core.KNOWN_DEC}, "
          f"site lat {core.KNOWN_SITE[0]}, lon {core.KNOWN_SITE[1]}, "
          f"height {core.KNOWN_SITE[2]:.0f} m (Paranal)")
    print(f"  computed: BJD_TDB {computed:.10f}")
    print(f"  expected: BJD_TDB {expected:.10f} "
          "(Astropy builtin ephemeris, independently reproduced by the "
          "geometric formula of Eastman et al. 2010 -- see bjd/core.py)")
    if diff_s <= core.KNOWN_TOLERANCE_S:
        print(f"  result  : PASS (difference {diff_s:.4f} s, "
              f"tolerance {core.KNOWN_TOLERANCE_S} s)")
        return 0
    print(f"  result  : FAIL (difference {diff_s:.4f} s exceeds "
          f"tolerance {core.KNOWN_TOLERANCE_S} s)")
    return 1


def _normalise_argv(argv):
    """Rewrite `--site -24.6,-70.4,2635` into `--site=-24.6,-70.4,2635`.

    Argparse mistakes a leading-minus coordinate string for a new flag
    and errors out with "expected one argument" -- which would break the
    documented usage for every southern-hemisphere observatory.
    """
    out = []
    i = 0
    while i < len(argv):
        if argv[i] == "--site" and i + 1 < len(argv):
            out.append("--site=" + argv[i + 1])
            i += 2
        else:
            out.append(argv[i])
            i += 1
    return out


def main(argv=None):
    parser = build_parser()
    raw = sys.argv[1:] if argv is None else list(argv)
    args = parser.parse_args(_normalise_argv(raw))
    if args.command is None:
        parser.print_help(sys.stderr)
        return 2
    try:
        if args.command == "convert":
            return cmd_convert(args)
        if args.command == "check":
            return cmd_check(args)
    except BjdError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # user errors must never show a traceback
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    sys.exit(main())
