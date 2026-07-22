from __future__ import annotations

import unittest

import pandas as pd

from external_validation import UnionFind, active_comparator_result, jader_age_20_49


class TestExternalValidation(unittest.TestCase):
    def test_union_find_transitive_links(self) -> None:
        uf = UnionFind()
        uf.union("3", "2")
        uf.union("2", "1")
        self.assertEqual(uf.find("3"), uf.find("1"))

    def test_jader_age_proxy(self) -> None:
        result = jader_age_20_49(pd.Series(["10歳代", "20歳代", "30歳代", "40歳代", "50歳代", "不明"]))
        self.assertEqual(result.tolist(), [False, True, True, True, False, False])

    def test_external_reconciliation(self) -> None:
        exposures = {
            "rivaroxaban": {"1", "2", "3"},
            "apixaban": {"3", "4", "5", "6"},
        }
        result = active_comparator_result(exposures, {"1", "4"})
        self.assertTrue(result["reconciliation_pass"])
        self.assertEqual(result["excluded_dual_exposure_reports"], 1)


if __name__ == "__main__":
    unittest.main()
