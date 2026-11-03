import tempfile
import unittest
import sys
from unittest.mock import patch
from pathlib import Path
import numpy as np
import pandas as pd
import annotate_sv326_with_encode4_ccre as a
import annotate_sv326_ccre_tissue_zscores as z
import annotate_sv326_ccre_target_genes as g


class AnnotationTests(unittest.TestCase):
    def test_single_base_and_invalid_coordinates(self):
        for typ in ("DEL", "INS", "VNTR"):
            row = pd.Series(a.parse_sv_id(f"chr1-10-10-{typ}-1"))
            self.assertEqual(a.sv_interval(row)[1:3], (9, 10))
        with self.assertRaises(ValueError):
            a.parse_sv_id("chr1-0-10-DEL-10")
        with self.assertRaises(ValueError):
            a.parse_sv_id("chr1-20-10-DEL-10")

    def test_exact_nearby_window(self):
        row = pd.Series(dict(SV_ID="chr1-10-10-INS-10", **a.parse_sv_id("chr1-10-10-INS-10")))
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "c.bed"
            p.write_text("chr1\t10\t11\tEH38D1\tEH38E1\tpELS\nchr1\t15\t16\tEH38D2\tEH38E2\tpELS\nchr1\t16\t17\tEH38D3\tEH38E3\tpELS\n")
            bed = a.read_ccre_bed(p)
            hits, _ = a.annotate_sv(row, {"chr1": bed}, 5, "span")
            self.assertEqual({h["distance_bp"] for h in hits}, {0, 5})
            self.assertTrue(all(h["relation"] == "nearby" for h in hits))
            hits, _ = a.annotate_sv(row, {"chr1": bed}, 0, "span")
            self.assertEqual(hits, [])

    def test_two_header_zscores(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "ENCFF560UZM.txt"
            p.write_text("rDHS\tENCSR000AAA\tENCSR000BBB\nrDHS\tprimary B cell\tliver\nEH38D1\t2.5\t0.2\n")
            rows, qc = z.scan_zscore_matrix(p, {"EH38D1"}, 1.64, 12, False)
            self.assertEqual(qc["header_rows"], 2)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["assay"], "H3K27ac")
            self.assertEqual(rows[0]["biosample"], "primary B cell")
            self.assertEqual(rows[0]["experiment_accession"], "ENCSR000AAA")
            self.assertEqual(rows[0]["zscore"], 2.5)
            rows, _ = z.scan_zscore_matrix(p, {"EH38D1"}, 1.64, 12, True)
            self.assertEqual(len(rows), 2)

    def test_join_accession_aliases_not_duplicates(self):
        hits = pd.DataFrame(dict(direct_hit_key=["x"], SV_ID=["sv"], ccre_id=["EH38E1"],
                                 e_accession=["EH38E1"], d_accession=["EH38D1"]))
        scores = pd.DataFrame(dict(ccre_accession=["EH38E1"], zscore=[2.0]))
        self.assertEqual(len(z.join_activity_to_hits(hits, scores)), 1)
        links = pd.DataFrame(dict(ccre_accession=["EH38E1"], gene_symbol=["GENE1"]))
        self.assertEqual(len(g.join_links_to_hits(hits, links)), 1)

    def test_bad_bed_not_silently_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "c.bed"
            p.write_text("chr1\t1.5\t10\tEH38D1\tEH38E1\tpELS\n")
            with self.assertRaises(ValueError):
                a.read_ccre_bed(p)

    def test_annotation_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pd.DataFrame(dict(SV_ID=["chr1-10-20-DEL-10"], PBS=[0.8])).to_csv(tmp/"sv.csv", index=False)
            (tmp/"ccre.bed").write_text("chr1\t9\t20\tEH38D1\tEH38E1\tpELS\n")
            with patch.object(sys, "argv", ["annotation", "--sv-csv", str(tmp/"sv.csv"), "--ccre-bed", str(tmp/"ccre.bed"), "--outdir", str(tmp/"ccre")]):
                a.main()
            hits = tmp/"ccre/sv326_encode4_ccre_hits.tsv"
            (tmp/"ENCFF560UZM.txt").write_text("rDHS\tENCSR000AAA\nrDHS\tprimary B cell\nEH38D1\t2.5\n")
            with patch.object(sys, "argv", ["zscore", "--ccre-hits", str(hits), "--zscore-files", str(tmp/"ENCFF560UZM.txt"), "--outdir", str(tmp/"z")]):
                z.main()
            link = tmp/"V4-hg38.Gene-Links.3D-Chromatin.txt"
            link.write_text("EH38E1\tENSG001\tGENE1\tprotein_coding\tHiC\tENCSR000AAA\tB cell\t1\t0.01\n")
            with patch.object(sys, "argv", ["links", "--ccre-hits", str(hits), "--gene-link-files", str(link), "--outdir", str(tmp/"g")]):
                g.main()
            result = pd.read_csv(tmp/"g/sv326_direct_ccre_gene_link_long.tsv", sep="\t")
            self.assertEqual(result.gene_symbol.tolist(), ["GENE1"])
            link.write_text("EH38E9\tENSG009\tGENE9\tprotein_coding\n")
            with patch.object(sys, "argv", ["links", "--ccre-hits", str(hits), "--gene-link-files", str(link), "--outdir", str(tmp/"empty")]):
                g.main()


if __name__ == "__main__":
    unittest.main()
