#!/usr/bin/env python3
"""
BEAST2 XML generation and execution for one 50-kb window.

Model actually used in the final analysis:
  HKY substitution, estimated base frequencies
  5 discrete gamma rate categories
  proportionInvariant FIXED at 0  (estimate="false", value="0.0")
  strict molecular clock
  calibrated Yule tree prior
  a single human-chimp calibration (LogNormal, meanInRealSpace=false)
  no gorilla
  1 chain x 50,000,000 states, log/tree every 1,000, storeEvery 5,000
  10% burn-in for summarisation

Run:
  python3 08_beast_analysis.py --xml  <write XML only>
  python3 08_beast_analysis.py --run  <write XML, run BEAST, TreeAnnotator>
"""
import argparse
import os
import subprocess
import sys

# ---- model constants (final production settings) ----
GAMMA_CATEGORIES = 5
PINV_VALUE = 0.0            # FIXED, not estimated
CHAIN_LENGTH = 50000000
LOG_EVERY = 1000
STORE_EVERY = 5000
BURNIN_FRAC = 0.10

XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast version="2.0" namespace="beast.base.evolution.alignment:beast.base.evolution.tree:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.branchratemodel:beast.base.evolution.speciation:beast.base.inference:beast.base.inference.operator:beast.base.inference.distribution:beast.base.math.distributions">
  <data id="data" spec="Alignment">
{sequences}
  </data>
  <run id="mcmc" spec="MCMC" chainLength="{chain_length}">
    <state id="state" storeEvery="{store_every}">
      <tree id="Tree.t" spec="beast.base.evolution.tree.Tree" name="stateNode">
        <taxonset id="TaxonSet" spec="TaxonSet">
          <alignment idref="data"/>
        </taxonset>
      </tree>
      <parameter id="kappa" name="stateNode" value="2.0" lower="0.0" estimate="true"/>
      <parameter id="pInv" name="stateNode" value="{pinv}" lower="0.0" upper="1.0" estimate="false"/>
      <parameter id="gammaShape" name="stateNode" value="0.1" lower="0.0" estimate="true"/>
      <parameter id="clockRate" name="stateNode" value="1.0e-9" lower="0.0" estimate="true"/>
    </state>

    <distribution id="posterior" spec="CompoundDistribution">
      <distribution id="treeLikelihood" spec="TreeLikelihood">
        <data idref="data"/>
        <tree idref="Tree.t"/>
        <siteModel id="SiteModel" spec="SiteModel" gammaCategoryCount="{gamma}">
          <parameter idref="pInv" name="proportionInvariant"/>
          <parameter idref="gammaShape" name="shape"/>
          <substModel id="hky" spec="HKY" kappa="@kappa"/>
        </siteModel>
        <branchRateModel id="StrictClock" spec="StrictClockModel" clock.rate="@clockRate"/>
      </distribution>
      <distribution id="CalibratedYule" spec="beast.base.evolution.speciation.CalibratedYuleModel"
                    tree="@Tree.t" birthRate="@birthRate"/>
    </distribution>

    <operator id="kappaScaler" spec="ScaleOperator" parameter="@kappa" scaleFactor="0.5" weight="0.1"/>
    <operator id="gammaScaler" spec="ScaleOperator" parameter="@gammaShape" scaleFactor="0.5" weight="0.1"/>
    <operator id="clockScaler" spec="ScaleOperator" parameter="@clockRate" scaleFactor="0.5" weight="0.1"/>
    <operator id="treeScaler" spec="ScaleOperator" parameter="@Tree.t" scaleFactor="0.5" weight="3.0"/>
    <operator id="treeRootScaler" spec="ScaleOperator" parameter="@Tree.t" scaleFactor="0.5" weight="3.0" rootOnly="true"/>

    <logger id="tracelog" fileName="{stem}.log" logEvery="{log_every}">
      <log idref="posterior"/>
      <log idref="treeLikelihood"/>
      <log idref="kappa"/>
      <log idref="gammaShape"/>
      <log idref="clockRate"/>
      <log id="mrcaHumanChimp" spec="beast.base.evolution.tree.MRCAPrior" tree="@Tree.t">
        <taxonset id="humanChimp" spec="TaxonSet">
          <taxon id="Chimpanzee" spec="Taxon"/>
          <taxon id="human" spec="Taxon"/>
        </taxonset>
      </log>
    </logger>

    <logger id="treelog" fileName="{stem}.trees" logEvery="{log_every}">
      <log idref="treeLikelihood"/>
    </logger>
  </run>
</beast>
"""


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_xml(fasta, out_xml):
    recs = read_fasta(fasta)
    seqs = "\n".join(
        f'    <sequence id="{n}" taxon="{n}" totalcount="{len(s)}" value="{s}"/>'
        for n, s in recs)
    stem = os.path.splitext(out_xml)[0]
    with open(out_xml, "w") as fh:
        fh.write(XML_TEMPLATE.format(
            sequences=seqs, chain_length=CHAIN_LENGTH,
            store_every=STORE_EVERY, pinv=PINV_VALUE,
            gamma=GAMMA_CATEGORIES, log_every=LOG_EVERY, stem=stem))
    return len(recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--xml", required=True)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--beast", default="beast")
    ap.add_argument("--treeannotator", default="treeannotator")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    n = write_xml(a.fasta, a.xml)
    print(f"wrote {a.xml} with {n} taxa")

    if not a.run:
        return

    stem = os.path.splitext(a.xml)[0]
    subprocess.run([a.beast, "-overwrite", "-beagle",
                    "-threads", str(a.threads),
                    "-seed", str(a.seed), a.xml], check=True)
    # 10% burn-in, median node heights
    subprocess.run([a.treeannotator, "-height", "median",
                    "-burnin", "10", "-limit", "0", "-topology", "MCC",
                    f"{stem}.trees", f"{stem}_MCC.tree"], check=True)


if __name__ == "__main__":
    main()
