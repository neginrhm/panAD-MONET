#!/usr/bin/env python3
"""Remove non-coding / provisional target genes from pySCENIC adjacency lists.

Applied to the sex-stratified GRNBoost2 outputs before `pyscenic ctx`
(roughly 42k of ~310k edges dropped per group in the original run).
A row is dropped when its `target` gene symbol matches one of the heuristic
patterns below (lincRNAs, antisense/intronic transcripts, clone-based
symbols, C#orf#, LOC#, pseudogene-like P# suffixes, microRNAs).

KNOWN ISSUE (kept as-is to reproduce the published regulons): the
pseudogene rule `P\d+$` also matches protein-coding genes whose symbol ends
in P + digits, e.g. SNAP25, AQP4, PLP1, SPP1, FOXP2, LAMP5, GRIP1, KCNIP4
and many ARHGAP* genes. In the dementia-male network this rule alone removed
19,545 of 302,347 edges (326 distinct target genes). Set
STRICT_PSEUDOGENE_RULE = True to restrict it to RPL/RPS-style pseudogenes
instead (this changes the results relative to the manuscript).

Usage:
    python 01_filter_adjacency_noncoding.py [files ...]
    (default: adata_{ct,mci,dm}_{male,female}_adj.csv in the working directory)

Output, per input <name>.csv:
    <name>_scenic_filteredGenes.csv
"""

import argparse
import glob
import re
from pathlib import Path

import pandas as pd

DEFAULT_INPUTS = [
    "adata_mci_male_adj.csv",
    "adata_mci_female_adj.csv",
    "adata_ct_female_adj.csv",
    "adata_ct_male_adj.csv",
    "adata_dm_male_adj.csv",
    "adata_dm_female_adj.csv",
]

COLUMN_WITH_GENE = "target"            # column tested for non-coding symbols
DROP_IF_EITHER_COLUMN_MATCHES = False  # True: also drop rows whose TF matches
STRICT_PSEUDOGENE_RULE = False         # False = original behavior (see KNOWN ISSUE above)

JUNK_PATTERNS = [
    r"^LINC",                        # long intergenic non-coding RNAs
    r"-AS\d*$", r"-AS$",             # antisense (e.g. EPHA1-AS1)
    r"-IT\d*$", r"-IT$",             # intronic transcripts
    r"^RP\d+-", r"^CTD-", r"^AC\d+", r"^AL\d+", r"^AP\d+", r"^BX\d+",  # clone-based symbols
    r"^C\d+orf\d+$",                 # C#orf# (provisional)
    r"^LOC\d+$",                     # LOC symbols
    (r"^RP[LS]\d+[A-Z]?\d*P\d+$" if STRICT_PSEUDOGENE_RULE
     else r"P\d+$"),                  # pseudogene-like P# suffix (e.g. RPL7P2)
    r"^MIR", r"^LET7",               # microRNAs
]
JUNK_RE = re.compile("|".join(JUNK_PATTERNS), flags=re.IGNORECASE)


def is_junk_gene(symbol) -> bool:
    if not isinstance(symbol, str) or symbol.strip() == "":
        return True
    return bool(JUNK_RE.search(symbol.strip()))


def filter_df(df: pd.DataFrame) -> pd.DataFrame:
    missing = {"TF", "target"} - set(df.columns)
    if missing:
        raise ValueError(f"input needs columns TF, target; missing: {missing}")
    if DROP_IF_EITHER_COLUMN_MATCHES:
        mask = df["TF"].map(is_junk_gene) | df["target"].map(is_junk_gene)
    else:
        mask = df[COLUMN_WITH_GENE].map(is_junk_gene)
    return df.loc[~mask].copy()


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="*", default=DEFAULT_INPUTS)
    args = p.parse_args()

    files = sorted({f for pat in args.files for f in glob.glob(pat)})
    if not files:
        raise SystemExit("no matching adjacency CSVs found")

    records = []
    for f in files:
        df = pd.read_csv(f)
        out = filter_df(df)
        out_path = str(Path(f).with_suffix("")) + "_scenic_filteredGenes.csv"
        out.to_csv(out_path, index=False)
        print(f"[OK] {Path(f).name}: {len(df)} -> {len(out)} kept ({len(df) - len(out)} dropped) -> {Path(out_path).name}")
        records.append({"file": Path(f).name, "rows_in": len(df), "rows_out": len(out),
                        "dropped": len(df) - len(out)})
    print(pd.DataFrame(records).sort_values("file").to_string(index=False))


if __name__ == "__main__":
    main()
