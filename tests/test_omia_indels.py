"""OMIA HGVS indels become VCF-style alleles with real letters from the species genome."""

import unittest

import inference_omia


class _Seq:
    def __init__(self, s):
        self.seq = s


class _Contig:
    def __init__(self, s):
        self.s = s

    def __getitem__(self, sl):
        return _Seq(self.s[sl])


GENOME = {"1": _Contig("ACGTACGTACGTACGT"), "2": _Contig("ACGTNNNNACGT")}  # chr1 1-based: 4 = T, 5 = A, 6 = C, 7 = G


def parse(hgvs):
    return inference_omia._extract_indel_from_hgvs(f"NC_000001.1:{hgvs}", GENOME, "1")


class OmiaIndelTest(unittest.TestCase):
    def test_deletion_is_anchored_on_the_preceding_letter(self):
        self.assertEqual(parse("g.5del"), (4, "TA", "T"))
        self.assertEqual(parse("g.5_7del"), (4, "TACG", "T"))

    def test_duplication_repeats_the_reference_letters(self):
        self.assertEqual(parse("g.5dup"), (5, "A", "AA"))

    def test_insertion_follows_the_anchor_letter(self):
        self.assertEqual(parse("g.5insTT"), (5, "A", "ATT"))

    def test_delins_replaces_the_span(self):
        self.assertEqual(parse("g.5_6delinsGG"), (5, "AC", "GG"))

    def test_unknown_inserted_sequence_is_skipped(self):
        for hgvs in ("g.5insN[226]", "g.5delinsN[226]", "g.5delAinsGG"):
            self.assertEqual(parse(hgvs), (None, None, None), hgvs)

    def test_declared_deleted_letters_must_match_the_genome(self):
        self.assertEqual(parse("g.5delA"), (4, "TA", "T"))
        self.assertEqual(parse("g.5delG"), (None, None, None))
        self.assertEqual(parse("g.5dupA"), (5, "A", "AA"))

    def test_insertion_range_must_be_adjacent(self):
        self.assertEqual(parse("g.5_6insTT"), (5, "A", "ATT"))
        self.assertEqual(parse("g.5_9insTT"), (None, None, None))

    def test_span_past_the_contig_end_is_skipped(self):
        self.assertEqual(parse("g.5_9000del"), (None, None, None))

    def test_assembly_gap_letters_are_skipped(self):
        got = inference_omia._extract_indel_from_hgvs("NC_000002.1:g.7del", GENOME, "2")
        self.assertEqual(got, (None, None, None))


if __name__ == "__main__":
    unittest.main()
