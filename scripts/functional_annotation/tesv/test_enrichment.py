import argparse
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import matched_enrichment as m
import prepare_regions as r


class EnrichmentTests(unittest.TestCase):
    def test_interval_boundary_nested_and_empty(self):
        sv = pd.DataFrame(dict(chrom=["chr1"]*4, start=[9, 10, 19, 20], end=[10, 11, 20, 21]))
        bed = pd.DataFrame(dict(chrom=["chr1", "chr1"], start=[10, 12], end=[20, 13]))
        self.assertEqual(m.overlap(sv, bed).tolist(), [False, True, True, False])
        self.assertFalse(m.overlap(sv, bed.iloc[:0]).any())

    def test_insertion_uses_reported_span_not_inserted_length(self):
        sv = m.parse_sv(["chr1-100-100-INS-500", "chr1-100-140-INS-40"], "span")
        self.assertEqual((sv.end-sv.start).tolist(), [1, 41])
        point = m.parse_sv(sv.SV_ID, "point")
        self.assertEqual((point.end-point.start).tolist(), [1, 1])

    def test_matching_sampling_no_leakage_reproducibility(self):
        sv = pd.DataFrame(dict(is_TESV=[True,True,False,False,False,False],
                               stratum=["a"]*6, length=[100]*6, MAF=[.2]*6,
                               interval_length=[100]*6, endpoint=[1,1,0,0,0,0]))
        with tempfile.TemporaryDirectory() as d:
            groups = m.matching(sv, Path(d))
        a, _, first = m.simulate(sv, ["endpoint"], groups, 20, 99)
        b, _, second = m.simulate(sv, ["endpoint"], groups, 20, 99)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(set(first)), 2)
        self.assertTrue((first >= 2).all())
        self.assertEqual(a.sum(), 0)

    def test_matching_shortage_fails_not_replacement(self):
        sv = pd.DataFrame(dict(is_TESV=[True, True, False], stratum=["a"]*3))
        with tempfile.TemporaryDirectory() as d, self.assertRaises(ValueError):
            m.matching(sv, Path(d))

    def test_ten_thousand_draws(self):
        sv = pd.DataFrame(dict(is_TESV=[True, True, False, False, False, False],
                               stratum=["a"]*6, length=[100]*6, MAF=[.2]*6,
                               interval_length=[100]*6, endpoint=[1,1,0,0,0,0]))
        with tempfile.TemporaryDirectory() as d:
            groups = m.matching(sv, Path(d))
        null, _, _ = m.simulate(sv, ["endpoint"], groups, 10000, 20260928)
        self.assertEqual(null.shape, (10000, 1))
        self.assertEqual((1 + (null[:, 0] >= 2).sum()) / 10001, 1 / 10001)

    def test_bh(self):
        np.testing.assert_allclose(m.bh([.01,.04,.03,1]), [.04,.053333333,.053333333,1])

    def test_reject_selected_only(self):
        with self.assertRaises(ValueError):
            r.reject_sv_conditioning(pd.DataFrame(columns=["interaction_id", "TESV"]))

    def test_study_exclusion_and_missing_maf(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            ids = ["chr1-10-20-DEL-10", "chr1-30-40-DEL-10"]
            pd.DataFrame(dict(svID_5caller=ids[:1]+["chr1-99-100-VNTR"], PBS=[.9,np.nan])).to_csv(d/"s.csv",index=False)
            u = pd.DataFrame(dict(KEY=ids, AF_TIB=[.1,.1], AF_HAN=[.3,.3]))
            u.to_csv(d/"u.csv",index=False)
            a = argparse.Namespace(selected=d/"s.csv",universe=d/"u.csv",ins_mode="span",expected_selected=1,
                                   maf_column="",tib_af_weight=.5,match_chrom=False,match_span=False)
            sv, _ = m.studies(a,d)
            self.assertEqual(sv.is_TESV.tolist(), [True,False])
            np.testing.assert_allclose(sv.MAF,[.2,.2])
            u.loc[1,"AF_HAN"] = np.nan
            u.to_csv(d/"u.csv",index=False)
            with self.assertRaises(ValueError):
                m.studies(a,d)

    def test_hic_complete_loop_not_other_loop(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            loops = pd.DataFrame(dict(interaction_id=["l1","l2"],chrom=["chr1"]*2,
                anchor1_start=[0,20000],anchor1_end=[20000,40000],anchor2_start=[60000,80000],anchor2_end=[80000,100000]))
            oe=loops.copy()
            for s in r.SAMPLES:
                oe[s]=[2 if s.startswith("A") else 0,0]
            oe.loc[0,"A-1-4"]=np.nan
            oe.to_csv(d/"oe.tsv",sep="\t",index=False)
            args=argparse.Namespace(hic_oe=d/"oe.tsv",hic_delta=1)
            r.hic_effect_regions(args,loops,d)
            self.assertEqual((d/"hic_effect.bed").stat().st_size,0)
            self.assertGreater((d/"hic_valid.bed").stat().st_size,0)

    def test_gene_links_use_all_genes_and_normalize_ids(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            rna=pd.DataFrame(dict(feature_id=["GENE1","ENSG0002"],gene_id=["GENE1","ENSG0002"],
                time6h=[1,0],time1d=[1,0],time3d=[1,0],time5d=[1,0],**{"adj.P.Val":[.01,.9]}))
            rna.to_csv(d/"rna.tsv",sep="\t",index=False)
            (d/"links.txt").write_text("EH38E111\tENSG0001.2\tGENE1\tprotein_coding\nEH38E222\tENSG0002.9\tGENE2\tprotein_coding\n")
            args=argparse.Namespace(rna_time=d/"rna.tsv",fdr=.05,min_logfc=.5,gene_link=[d/"links.txt"])
            self.assertEqual(r.responsive_accessions(args,d),{"EH38E111"})

    def test_blocked_sparse_hic_pixel_calculation(self):
        class Block:
            def __init__(self,x): self.x=x
            def tocsr(self): return self
            def __getitem__(self,key): return self.x[key]
        class Matrix:
            def __getitem__(self,key):
                a=np.ones((5,5))*4
                a[1,4]=0
                return Block(a[key])
        class Cooler:
            binsize=20000
            chromnames=["chr1"]
            def __init__(self,path): pass
            def extent(self,chrom): return 0,5
            def bins(self): return self
            def fetch(self,chrom): return pd.DataFrame(dict(weight=[1]*5))
            def matrix(self,**kwargs): return Matrix()
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            (d/"x.mcool").write_text("mock")
            pd.DataFrame(dict(region1=["chr1"],dist=[3],**{"balanced.avg.smoothed":[2]})).to_csv(d/"e.tsv",sep="\t",index=False)
            pd.DataFrame(dict(sample=r.SAMPLES,mcool=[str(d/"x.mcool")]*8,expected_cis=[str(d/"e.tsv")]*8)).to_csv(d/"meta.tsv",sep="\t",index=False)
            loops=pd.DataFrame(dict(interaction_id=["a","b"],chrom=["chr1"]*2,
                anchor1_start=[0,20000],anchor1_end=[20000,40000],anchor2_start=[60000,80000],anchor2_end=[80000,100000]))
            with patch.dict(sys.modules,{"cooler":types.SimpleNamespace(Cooler=Cooler)}):
                result=r.quantify_hic(loops,d/"meta.tsv",d,2)
            self.assertEqual(result.loc[0,"A-1-1"],1)
            self.assertTrue(pd.isna(result.loc[1,"A-1-1"]))

    def test_all_endpoints_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            peak = pd.DataFrame(dict(peak_id=["p1", "p2"], chr=["chr1"]*2, start=[100,40000], end=[150,40050]))
            for pop in ["A", "W"]:
                for donor in range(1,6):
                    for t in range(1,6): peak[f"{pop}-{donor}-{t}"] = [3,4]
            peak.to_csv(d/"rpkm.tsv",sep="\t",index=False)
            pd.DataFrame(dict(feature_id=["p1","p2"],logFC=[1,0],**{"adj.P.Val":[.01,.9]})).to_csv(d/"dor.tsv",sep="\t",index=False)
            loops=pd.DataFrame(dict(interaction_id=["l1","l2"],chrom=["chr1"]*2,
                anchor1_start=[0,40000],anchor1_end=[20000,60000],anchor2_start=[60000,100000],anchor2_end=[80000,120000]))
            loops.to_csv(d/"loops.tsv",sep="\t",index=False)
            oe=loops.copy()
            for s in r.SAMPLES: oe[s] = [2 if s.startswith("A") else 0,0]
            oe.to_csv(d/"oe.tsv",sep="\t",index=False)
            pd.DataFrame(dict(feature_id=["G1"],gene_id=["G1"],time6h=[1],time1d=[1],time3d=[1],time5d=[1],**{"adj.P.Val":[.01]})).to_csv(d/"rna.tsv",sep="\t",index=False)
            (d/"links.txt").write_text("EH38E111\tENSG001\tG1\tprotein_coding\n")
            (d/"ccre.bed").write_text("chr1\t100\t150\tEH38D111\tEH38E111\tpELS\n")
            argv=["prepare_regions.py","--atac-rpkm",str(d/"rpkm.tsv"),"--atac-population",str(d/"dor.tsv"),
                "--hic-loops",str(d/"loops.tsv"),"--hic-oe",str(d/"oe.tsv"),"--rna-time",str(d/"rna.tsv"),
                "--ccre-bed",str(d/"ccre.bed"),"--gene-link",str(d/"links.txt"),"--outdir",str(d/"regions")]
            with patch.object(sys,"argv",argv): r.main()
            ids=[f"chr1-{i}-{i+20}-DEL-20" for i in [110,40010,2000,3000,80000,90000,100010,500000]]
            pd.DataFrame(dict(KEY=ids,AF_TIB=[.2]*8,AF_HAN=[.2]*8)).to_csv(d/"universe.csv",index=False)
            pd.DataFrame(dict(svID_5caller=ids[:2],PBS=[.8,.9])).to_csv(d/"selected.csv",index=False)
            argv=["matched_enrichment.py","--selected",str(d/"selected.csv"),"--universe",str(d/"universe.csv"),
                "--manifest",str(d/"regions/region_manifest.tsv"),"--outdir",str(d/"results"),
                "--expected-selected","2","--n-perm","19","--require-complete"]
            with patch.object(sys,"argv",argv): m.main()
            result=pd.read_csv(d/"results/enrichment_summary.tsv",sep="\t").set_index("endpoint")
            self.assertTrue(result.status.eq("READY").all())
            self.assertEqual(len(result),9)
            self.assertEqual(result.loc["ATAC_open_direct","observed_SVs"],2)
            self.assertEqual(result.loc["ATAC_population_DOR_direct","observed_SVs"],1)
            self.assertEqual(result.loc["HiC_population_effect_direct","observed_SVs"],1)
            self.assertEqual(result.loc["cCRE_hypoxia_gene_direct","observed_SVs"],1)
            self.assertTrue(result.FDR_enrichment.between(0,1).all())


if __name__ == "__main__":
    unittest.main()
