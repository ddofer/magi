"""Rule layer: provenance columns and the exported rule table."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "clinvar_pipeline/vendor"))
from nt_mechanism_assignment import assign_mechanisms_vectorized, get_hypothesis_logic, rule_table


class RuleLayerTest(unittest.TestCase):
    def test_provenance_columns(self):
        df = pd.DataFrame({
            "D_BED_exon": [0.2, np.nan, 0.01],
            "D_BED_intron": [-0.12, np.nan, 0.0],
            "LLR": [0.0, np.nan, -0.1],
            "MLM_Prior": [0.5, np.nan, 0.5],
        })
        out = assign_mechanisms_vectorized(df, get_hypothesis_logic("snp"))
        self.assertEqual(list(out["NT_Mechanism"]),
                         ["SPLICE_INDUCED_EXONIZATION", "UNCERTAIN_SIGNIFICANCE", "UNCERTAIN_SIGNIFICANCE"])
        self.assertEqual(list(out["NT_Tier"]), [1, -1, -1])
        self.assertEqual(out["NT_Rule"][0], "exon>=0.15 [obs 0.200] AND intron<=-0.10 [obs -0.120]")
        self.assertEqual(list(out["NT_Scored"]), [True, False, True])

    def test_only_fired_primaries_listed(self):
        df = pd.DataFrame({"D_BED_splice_donor": [-0.3, -0.3], "D_BED_splice_acceptor": [0.0, -0.25]})
        out = assign_mechanisms_vectorized(df, get_hypothesis_logic("snp"))
        self.assertEqual(list(out["NT_Rule"]), [
            "splice_donor<=-0.20 [obs -0.300]",
            "splice_donor<=-0.20 [obs -0.300] OR splice_acceptor<=-0.20 [obs -0.250]",
        ])
        self.assertEqual(list(out["NT_Tier"]), [0, 0])

    def test_rule_table(self):
        for vt in ("snp", "indel"):
            table = rule_table(vt).drop_duplicates("mechanism").set_index("mechanism")
            self.assertEqual(len(table), len(get_hypothesis_logic(vt)))
            self.assertEqual(list(table.index[table["can_be_assigned"] == False]),
                             ["RECIPROCAL_EXON_GAIN_INTRON_LOSS"])
            self.assertEqual(table.loc[["SPLICE_SITE_DESTROYED", "START_LOSS", "STOP_CODON_DISRUPTION"],
                                       "order_within_tier"].tolist(), [0, 1, 2])


if __name__ == "__main__":
    unittest.main()
