from __future__ import annotations

import unittest

from pv_stats import benjamini_hochberg, log_or_heterogeneity, minimum_detectable_or, two_by_two_statistics
from faers_full import fda_date_to_period


class TestPvStats(unittest.TestCase):
    def test_effect_measures(self) -> None:
        result = two_by_two_statistics(20, 80, 5, 95)
        self.assertAlmostEqual(result["ror"], 4.75)
        self.assertAlmostEqual(result["prr"], 4.0)
        self.assertGreater(result["ror_ci95_upper"], result["ror"])
        self.assertLess(result["ror_ci95_lower"], result["ror"])

    def test_exact_interval_handles_zero_cell(self) -> None:
        result = two_by_two_statistics(2, 10, 0, 12)
        self.assertTrue(result["jeffreys_0_5_correction"])
        self.assertGreaterEqual(result["conditional_exact_ci95_lower"], 0)
        self.assertEqual(result["conditional_exact_ci95_upper"], float("inf"))

    def test_log_or_heterogeneity(self) -> None:
        result = log_or_heterogeneity(
            [("hcp", 88, 1188, 27, 618), ("consumer", 481, 1928, 39, 1020)]
        )
        self.assertEqual(result["degrees_freedom"], 1)
        self.assertLess(result["p_heterogeneity"], 0.05)

    def test_bh_is_monotone_after_ranking(self) -> None:
        p = [0.04, 0.001, 0.03, 0.20]
        q = benjamini_hochberg(p)
        ordered = sorted(zip(p, q))
        self.assertTrue(all(ordered[i][1] <= ordered[i + 1][1] for i in range(len(ordered) - 1)))
        self.assertTrue(all(0 <= value <= 1 for value in q))

    def test_mde_is_finite(self) -> None:
        value = minimum_detectable_or(500, 500, 25)
        self.assertGreater(value, 1.0)
        self.assertLess(value, 10.0)

    def test_fda_date_quarter_assignment(self) -> None:
        self.assertEqual(fda_date_to_period(20210430), "2021Q2")
        self.assertEqual(fda_date_to_period("20151201"), "2015Q4")
        with self.assertRaises(ValueError):
            fda_date_to_period(20211301)


if __name__ == "__main__":
    unittest.main()
