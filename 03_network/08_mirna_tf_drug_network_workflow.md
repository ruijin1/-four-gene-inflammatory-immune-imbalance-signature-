# 08  mRNA-miRNA / mRNA-TF / mRNA-drug regulatory networks

Manuscript: Methods, "Networks ..."; Fig. 9A-C

## 8.1 mRNA-miRNA (miRTarBase v9.0)
- Query: https://mirtarbase.cuhk.edu.cn/ , species hsa, input the 4 hub genes
- Keep only records with **experimental support** (strong / weak evidence)
- Export `mirtarbase_hub_pairs.csv` (miRNA - target - support type)
- Build the network in Cytoscape 3.10.3 (Fig. 9A)

## 8.2 mRNA-TF (JASPAR v2024)
- Query: https://jaspar.genereg.net/ , profile = vertebrates (Homo sapiens)
- Take the predicted binding sites within the promoter-enhancer region of each hub gene
  (default +/- 2 kb)
- Relative score threshold: default 80%; export `jaspar_TF_pairs.csv`
- Build the network in Cytoscape (Fig. 9B)

## 8.3 mRNA-drug (DGIdb v5.0)
- Query: https://www.dgidb.org/ , input the 4 hub genes
- Include both **direct and indirect** targeting relationships
- Export `dgidb_result.tsv` (gene - drug - interaction type - source)
- Final result: **44 drug-target pairs**; build the network in Cytoscape (Fig. 9C)

## Note on reproducibility

All three databases offer **bulk downloads**: the raw export files (csv / tsv) can be uploaded to
the code repository together with the Cytoscape session. GUI operations cannot be reproduced from a
script alone, so **the raw export files are the key to reproducibility**.
