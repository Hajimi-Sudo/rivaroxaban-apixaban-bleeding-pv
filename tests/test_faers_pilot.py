from __future__ import annotations

import math
import json
import unittest

import pandas as pd

from faers_pilot import TABLE_PATTERNS, age_in_years, deduplicate_demo_frame, reporting_odds_ratio


class TestFaersPilot(unittest.TestCase):
    def test_age_conversion(self) -> None:
        age = pd.Series([24, 240, 3, 365, 4, 10])
        code = pd.Series(["YR", "MON", "DEC", "DY", "WK", "BAD"])
        converted = age_in_years(age, code)
        self.assertAlmostEqual(converted.iloc[0], 24.0)
        self.assertAlmostEqual(converted.iloc[1], 20.0)
        self.assertAlmostEqual(converted.iloc[2], 30.0)
        self.assertAlmostEqual(converted.iloc[3], 365 / 365.25)
        self.assertAlmostEqual(converted.iloc[4], 28 / 365.25)
        self.assertTrue(pd.isna(converted.iloc[5]))

    def test_ror_formula(self) -> None:
        result = reporting_odds_ratio(20, 80, 5, 95)
        self.assertAlmostEqual(result["ror"], 4.75)
        self.assertLess(result["ci95_lower"], result["ror"])
        self.assertGreater(result["ci95_upper"], result["ror"])
        self.assertFalse(result["jeffreys_0_5_correction"])

    def test_zero_cell_correction(self) -> None:
        result = reporting_odds_ratio(0, 10, 2, 20)
        self.assertTrue(math.isfinite(result["ror"]))
        self.assertTrue(result["jeffreys_0_5_correction"])
        json.dumps(result)

    def test_2018q1_new_demo_filename(self) -> None:
        self.assertIsNotNone(TABLE_PATTERNS["DEMO"].search("ascii/DEMO18Q1_new.txt"))
        self.assertIsNone(TABLE_PATTERNS["DEMO"].search("ascii/demo18q1.pdf"))

    def test_same_primaryid_corrected_fda_date_keeps_latest(self) -> None:
        frame = pd.DataFrame(
            [
                {"primaryid": "1001", "caseid": "100", "caseversion": "1", "fda_dt": "20211102", "age": "30", "age_cod": "YR", "sex": "F", "period": "2021Q4"},
                {"primaryid": "1001", "caseid": "100", "caseversion": "1", "fda_dt": "20220301", "age": "30", "age_cod": "YR", "sex": "F", "period": "2022Q1"},
            ]
        )
        result = deduplicate_demo_frame(frame)
        self.assertEqual(len(result), 1)
        self.assertEqual(int(result.iloc[0]["fda_dt_num"]), 20220301)
        self.assertEqual(len(result.attrs["primaryid_fda_dt_conflicts"]), 2)

    def test_same_primaryid_conflicting_caseversion_rejected(self) -> None:
        frame = pd.DataFrame(
            [
                {"primaryid": "1001", "caseid": "100", "caseversion": "1", "fda_dt": "20211102", "age": "30", "age_cod": "YR", "sex": "F", "period": "2021Q4"},
                {"primaryid": "1001", "caseid": "100", "caseversion": "2", "fda_dt": "20220301", "age": "30", "age_cod": "YR", "sex": "F", "period": "2022Q1"},
            ]
        )
        with self.assertRaises(ValueError):
            deduplicate_demo_frame(frame)


if __name__ == "__main__":
    unittest.main()
