# 07  PPI network and hub-gene selection (STRING + Cytoscape/CytoHubba)

Manuscript: Methods, "Networks of protein-protein interactions"; Fig. 4A-B

> This step is a **Cytoscape graphical-interface procedure**; no script can be archived. The
> complete parameters and steps needed to reproduce it are recorded below.

## Steps

1. **Input**: gene symbols of the 57 DE-IRGs (`results/DE_IRG_57.csv`).
2. **STRING v12.0** (https://string-db.org)
   - species *Homo sapiens*; `minimum interaction score = 0.400` (moderate-to-high confidence)
   - export TSV (node1 / node2 / combined_score) -> **48 connected genes / 172 edges**
   - sort by `combined_score` descending; also export `string_interactions.tsv` for audit
3. Import the edge table into **Cytoscape v3.10.3** -> PPI network (Fig. 4A).
4. **CytoHubba plugin**: take the **Top 30** genes by each of the **EPC**, **MCC** and **Degree**
   algorithms.
5. Intersect the three Top-30 lists -> **4 core hub genes**: CD3G, PRTN3, CD6, THBS1
   (Fig. 4B, Venn diagram).
6. **Export**: `CytoHubba_top30_EPC.txt` / `_MCC.txt` / `_Degree.txt`, and the intersected gene
   list.

## Files released alongside this workflow

- `string_interactions.tsv` (raw edge table with combined_score)
- `CytoHubba_top30_*.txt` (one file per algorithm)
- `cytoscape_session.cys` (optional; lets reviewers reproduce the network layout)
