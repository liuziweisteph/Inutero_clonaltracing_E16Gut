# Run from project root. Requires E9E10/E9E10NC.AC_neural_and_epi.h5ad
# (download or generate via E9E10/E9E10_analysis).

import re
import scanpy as sc
import pandas as pd
import numpy as np

INPUT_FILE = "E9E10/E9E10NC.AC_neural_and_epi.h5ad"
OUTPUT_FILE = "Automatic_annotation/scanpy_E9E10NC.AC_harmony_clusters.csv"

NC_AC = sc.read_h5ad(INPUT_FILE)
sc.tl.rank_genes_groups(NC_AC, groupby="leiden_res_0.40", method="wilcoxon", use_raw=False)
markers_res = NC_AC.uns["rank_genes_groups"]


def sort_groups(groups):
    """Sort cluster/group names numeric-first (matches SCSA/Seurat ordering)."""
    def key(g):
        s = str(g)
        return (0, int(s)) if s.isdigit() else (1, s)
    return sorted(groups, key=key)


def format_scsa_results(result_dict, output_file):
    """
    Format scanpy differential expression results into SCSA CSV format,
    keeping only positive markers (logFC > 1, adj p < 0.05).
    """
    names = result_dict["names"]
    scores = result_dict["scores"]
    pvals_adj = result_dict["pvals_adj"]
    logfoldchanges = result_dict["logfoldchanges"]

    groups = sort_groups(names.dtype.names)
    n_clusters = len(groups)

    # Build header per cluster: _n (gene), _l (logFC), _s (score), _p (adj p)
    header_cols = []
    for group in groups:
        header_cols.extend([
            f"{group}_n",
            f"{group}_l",
            f"{group}_s",
            f"{group}_p",
        ])

    positive_indices = {}
    max_len = 0
    for group in groups:
        lfc = logfoldchanges[group]
        padj = pvals_adj[group]
        # Boolean mask
        mask = (padj < 0.05) & (lfc > 1)
        idx = np.nonzero(mask)[0]
        positive_indices[group] = idx
        max_len = max(max_len, len(idx))

    data_rows = []
    for rank in range(max_len):
        row = [rank]  # first column is rank index
        for group in groups:
            idxs = positive_indices[group]
            if rank < len(idxs):
                i = idxs[rank]
                gene_name = names[group][i]
                gene_name_clean = re.sub(r"[^\x00-\x7F]+", "", str(gene_name)).encode("ascii", "ignore").decode("ascii")
                lfc = logfoldchanges[group][i]
                score = scores[group][i]
                padj = pvals_adj[group][i]
                row.extend([gene_name_clean, lfc, score, padj])
            else:
                row.extend(["", np.nan, np.nan, np.nan])
        data_rows.append(row)

    columns = [""] + header_cols
    df = pd.DataFrame(data_rows, columns=columns)
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    print(f"- Groups: {groups}")
    print(f"- Positive markers max rows: {max_len}")
    return df


if __name__ == "__main__":
    format_scsa_results(markers_res, OUTPUT_FILE)
    print("Marker analysis complete. Results saved to", OUTPUT_FILE)
