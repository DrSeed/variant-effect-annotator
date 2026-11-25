# Self-contained demo: every possible SNV in a simulated gene, plus a small "VCF" of mixed variants.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from annotate import CODON, annotate, annotate_snv

SENSE = [c for c, aa in CODON.items() if aa != "*" and c != "ATG"]


def random_cds(n_codons, rng):
    return "ATG" + "".join(rng.choice(SENSE, n_codons - 2)) + "TAA"


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(23)
    cds = random_cds(500, rng)
    every = pd.DataFrame([annotate_snv(cds, p, alt) for p in range(len(cds)) for alt in "ACGT" if alt != cds[p]])

    vcf = []
    for _ in range(25):
        p = int(rng.integers(3, len(cds) - 6))
        kind = rng.random()
        if kind < 0.7:
            vcf.append((p, cds[p], rng.choice([b for b in "ACGT" if b != cds[p]])))
        elif kind < 0.85:
            vcf.append((p, cds[p], cds[p] + "".join(rng.choice(list("ACGT"), int(rng.integers(1, 4))))))
        else:
            k = int(rng.integers(1, 4))
            vcf.append((p, cds[p:p + k + 1], cds[p]))
    calls = pd.DataFrame(annotate(cds, sorted(vcf)))
    calls.to_csv("results/annotated_variants.csv", index=False)

    table = every.groupby(["codon_position", "effect"]).size().unstack(fill_value=0)
    table = (100 * table.div(table.sum(axis=1), axis=0)).round(1)
    table.to_csv("results/summary.csv")

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    colours = {"synonymous": "#2a7f62", "missense": "#ff7f0e", "stop_gained": "#c0504d", "stop_lost": "#9467bd", "start_lost": "#444444"}
    bottom = np.zeros(3)
    for eff, c in colours.items():
        if eff in table:
            axes[0].bar(["1st", "2nd", "3rd"], table[eff], bottom=bottom, color=c, label=eff)
            bottom += table[eff].values
    axes[0].set_xlabel("codon position of the SNV")
    axes[0].set_ylabel("% of all possible SNVs")
    axes[0].set_title("Third positions are mostly silent")
    axes[0].legend(fontsize=8)

    third = every[every["codon_position"] == 3]
    tt = (100 * pd.crosstab(third["type"], third["effect"] == "synonymous", normalize="index")).round(1)
    axes[1].bar(tt.index, tt[True], color="#2a7f62")
    axes[1].set_ylabel("% synonymous at 3rd position")
    axes[1].set_title("Transitions are silent more often")

    miss = every[every["effect"] == "missense"]
    cls = (miss["note"].where(miss["note"] == "conservative", "class-changing")).value_counts()
    axes[2].bar(cls.index, cls.values, color=["#ff7f0e", "#bbbbbb"][:len(cls)])
    axes[2].set_ylabel("possible missense SNVs")
    axes[2].set_title("Most missense changes alter residue class")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table)
    print(calls[["c_change", "effect", "p_change", "note"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
