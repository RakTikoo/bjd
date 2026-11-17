"""Tests for bjd.

Run with the project venv:  ../.venv/bin/python run_tests.py
(or any Python with astropy installed, from the repo root.)
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from astropy.table import Table  # noqa: E402

from bjd import core  # noqa: E402

PARANAL = "-24.6272,-70.4042,2635"
RA, DEC = core.KNOWN_RA, core.KNOWN_DEC


def run_cli(*argv):
    """Run the CLI in a subprocess exactly as a user would."""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "bjd", *argv],
        capture_output=True, text=True, cwd=ROOT, env=env)


class KnownValueTests(unittest.TestCase):
    def test_known_value(self):
        # Reference case (see bjd/core.py for the full citation):
        #   JD 2460000.5 UTC, RA 84.291188, Dec -20.273999,
        #   Paranal (lat -24.6272, lon -70.4042, 2635 m).
        # Expected BJD_TDB = 2460000.5019914582, computed with Astropy
        # 8.0.2 (builtin ephemeris) and independently reproduced to
        # <1 ms by the direct geometric formula LTT = (r_obs . n_hat)/c
        # of Eastman, Siverd & Gaudi (2010, PASP, 122, 935).
        computed, expected = core.run_known_case()
        diff_s = abs(computed - expected) * 86400.0
        self.assertLess(diff_s, 0.05,
                        f"BJD_TDB off by {diff_s:.4f} s")

    def test_correction_magnitude_plausible(self):
        # Round-trip sanity: the total shift from a JD_UTC reading to
        # BJD_TDB combines TDB-UTC (~69 s) with the light-travel time
        # (at most ~500 s, one AU at light speed). For the reference
        # case it must land in the plausible tens-to-hundreds range.
        location = core.parse_site(PARANAL)
        bjd, times = core.to_bjd_tdb([core.KNOWN_INPUT_JD], "jd", "utc",
                                     RA, DEC, location)
        shift_s = float((bjd.jd[0] - times.jd[0]) * 86400.0)
        self.assertGreater(shift_s, 30.0)
        self.assertLess(shift_s, 600.0)
        ltt_s = shift_s - float((times.tdb.jd[0] - times.jd[0]) * 86400.0)
        self.assertGreater(abs(ltt_s), 1.0)
        self.assertLess(abs(ltt_s), 500.0)


class FormatTests(unittest.TestCase):
    def setUp(self):
        self.location = core.parse_site(PARANAL)

    def test_mjd_input(self):
        # MJD 60000.0 is the same instant as JD 2460000.5.
        bjd, _ = core.to_bjd_tdb([60000.0], "mjd", "utc",
                                 RA, DEC, self.location)
        diff_s = abs(float(bjd.jd[0]) - core.KNOWN_BJD_TDB) * 86400.0
        self.assertLess(diff_s, 0.05)

    def test_iso_input(self):
        # 2023-02-25 00:00:00 UTC is the same instant as JD 2460000.5.
        bjd, _ = core.to_bjd_tdb(["2023-02-25 00:00:00"], "iso", "utc",
                                 RA, DEC, self.location)
        diff_s = abs(float(bjd.jd[0]) - core.KNOWN_BJD_TDB) * 86400.0
        self.assertLess(diff_s, 0.05)

    def test_tt_scale_input(self):
        # The same instant stated in TT: JD_TT = JD_UTC + 69.184 s.
        jd_tt = core.KNOWN_INPUT_JD + 69.184 / 86400.0
        bjd, _ = core.to_bjd_tdb([jd_tt], "jd", "tt",
                                 RA, DEC, self.location)
        diff_s = abs(float(bjd.jd[0]) - core.KNOWN_BJD_TDB) * 86400.0
        self.assertLess(diff_s, 0.05)


class RefusalTests(unittest.TestCase):
    """The heart of the tool: no --format or --scale, no conversion."""

    def _convert(self, tmp, *extra):
        csv_path = Path(tmp) / "t.csv"
        csv_path.write_text("jd_utc\n2460000.5\n")
        return run_cli("convert", str(csv_path), "--time-col", "jd_utc",
                       "--ra", str(RA), "--dec", str(DEC),
                       "--site", PARANAL, *extra)

    def test_missing_scale_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = self._convert(tmp, "--format", "jd")
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("refuses to guess", proc.stderr)
        self.assertIn("--scale", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_missing_format_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = self._convert(tmp, "--scale", "utc")
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("refuses to guess", proc.stderr)
        self.assertIn("--format", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)


class EndToEndTests(unittest.TestCase):
    def test_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "times.csv"
            src.write_text("jd_utc,note\n2460000.5,a\n2460000.51,b\n")
            out = Path(tmp) / "out.csv"
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL, "--out", str(out))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("median shift applied", proc.stdout)
            table = Table.read(out)
            self.assertIn("bjd_tdb", table.colnames)
            self.assertEqual(len(table), 2)
            self.assertEqual(list(table["note"]), ["a", "b"])
            diff_s = abs(float(table["bjd_tdb"][0])
                         - core.KNOWN_BJD_TDB) * 86400.0
            self.assertLess(diff_s, 0.05)

    def test_fits_input_to_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "times.fits"
            Table({"jd_utc": [2460000.5, 2460000.51]}).write(src)
            out = Path(tmp) / "out.csv"
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL, "--out", str(out))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            table = Table.read(out)
            self.assertIn("bjd_tdb", table.colnames)
            self.assertEqual(len(table), 2)

    def test_check_subcommand(self):
        proc = run_cli("check")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("PASS", proc.stdout)


class AdversarialTests(unittest.TestCase):
    """Bad inputs must produce clean one-line errors, never tracebacks."""

    def _expect_clean_error(self, proc):
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("error", proc.stderr.lower())
        self.assertNotIn("Traceback", proc.stderr)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "empty.csv"
            src.write_text("")
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL)
        self._expect_clean_error(proc)

    def test_missing_time_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "t.csv"
            src.write_text("mjd_utc\n60000.0\n")
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL)
        self._expect_clean_error(proc)
        self.assertIn("not found", proc.stderr)

    def test_malformed_site(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "t.csv"
            src.write_text("jd_utc\n2460000.5\n")
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", "abc,def")
        self._expect_clean_error(proc)
        self.assertIn("malformed", proc.stderr)

    def test_named_site_explained(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "t.csv"
            src.write_text("jd_utc\n2460000.5\n")
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", "paranal")
        self._expect_clean_error(proc)
        self.assertIn("offline", proc.stderr)

    def test_site_lat_lon_swapped_is_caught_when_impossible(self):
        # lat=-120 cannot be a latitude: clean error, not a silent swap.
        with self.assertRaises(core.BjdError):
            core.parse_site("-120.0,-70.4,2635")


class ReviewRegressionTests(unittest.TestCase):
    """Regressions from the 2026-11-17 adversarial review (see REVIEW.md)."""

    def test_site_height_nan_refused(self):
        # Previously: rc=0 and an all-nan bjd_tdb column -- silent corruption.
        with self.assertRaises(core.BjdError):
            core.parse_site("-24.6272,-70.4042,nan")

    def test_site_height_inf_refused(self):
        with self.assertRaises(core.BjdError):
            core.parse_site("-24.6272,-70.4042,inf")

    def test_site_height_implausible_refused(self):
        with self.assertRaises(core.BjdError):
            core.parse_site("-24.6272,-70.4042,1e9")
        with self.assertRaises(core.BjdError):
            core.parse_site("-24.6272,-70.4042,-100000")

    def test_site_height_plausible_edges_accepted(self):
        # Dead Sea shore, highest observatories, stratospheric balloon.
        for h in ("-430", "5640", "40000"):
            loc = core.parse_site(f"-24.6272,-70.4042,{h}")
            self.assertAlmostEqual(float(loc.height.to_value("m")),
                                   float(h))

    def test_iso_t_separator_accepted(self):
        location = core.parse_site(PARANAL)
        bjd, _ = core.to_bjd_tdb(["2023-02-25T00:00:00"], "iso", "utc",
                                 RA, DEC, location)
        diff_s = abs(float(bjd.jd[0]) - core.KNOWN_BJD_TDB) * 86400.0
        self.assertLess(diff_s, 0.05)

    def test_iso_trailing_z_refused_cleanly(self):
        location = core.parse_site(PARANAL)
        with self.assertRaises(core.BjdError):
            core.to_bjd_tdb(["2023-02-25T00:00:00Z"], "iso", "utc",
                            RA, DEC, location)

    def test_api_bad_format_and_scale_are_bjd_errors(self):
        location = core.parse_site(PARANAL)
        with self.assertRaises(core.BjdError):
            core.to_bjd_tdb([2460000.5], "bogus", "utc",
                            RA, DEC, location)
        with self.assertRaises(core.BjdError):
            core.to_bjd_tdb([2460000.5], "jd", "bogus",
                            RA, DEC, location)

    def test_api_bad_target_is_bjd_error_not_nan(self):
        location = core.parse_site(PARANAL)
        with self.assertRaises(core.BjdError):
            core.to_bjd_tdb([2460000.5], "jd", "utc",
                            RA, 120.0, location)
        with self.assertRaises(core.BjdError):
            core.to_bjd_tdb([2460000.5], "jd", "utc",
                            float("nan"), DEC, location)

    def test_bom_csv_converts(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "bom.csv"
            src.write_bytes(b"\xef\xbb\xbfjd_utc,note\n2460000.5,a\n")
            out = Path(tmp) / "out.csv"
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL, "--out", str(out))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("bjd_tdb", Table.read(out).colnames)

    def test_out_same_as_input_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "t.csv"
            src.write_text("jd_utc\n2460000.5\n")
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL, "--out", str(src))
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("input file itself", proc.stderr)
            self.assertEqual(src.read_text(), "jd_utc\n2460000.5\n")

    def test_corrupt_file_error_is_one_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "bad.fits"
            src.write_bytes(b"NOT A FITS FILE AT ALL " * 10)
            proc = run_cli("convert", str(src), "--time-col", "jd_utc",
                           "--format", "jd", "--scale", "utc",
                           "--ra", str(RA), "--dec", str(DEC),
                           "--site", PARANAL)
            self.assertNotEqual(proc.returncode, 0)
            self.assertEqual(len(proc.stderr.strip().splitlines()), 1)
            self.assertNotIn("Traceback", proc.stderr)


if __name__ == "__main__":
    unittest.main()
