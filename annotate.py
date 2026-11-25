# Minimal coding-variant annotator: SNVs and indels in a CDS, HGVS-style protein notation.
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly",
         "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
         "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}
CLASS = {**dict.fromkeys("AVLIMFWC", "hydrophobic"), **dict.fromkeys("STNQGPY", "polar/small"),
         **dict.fromkeys("DE", "acidic"), **dict.fromkeys("KRH", "basic")}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def annotate_snv(cds, pos, alt):
    """pos is 0-based within the CDS. Returns a dict describing the consequence."""
    ref = cds[pos]
    codon_i, offset = divmod(pos, 3)
    ref_codon = cds[codon_i * 3: codon_i * 3 + 3]
    alt_codon = ref_codon[:offset] + alt + ref_codon[offset + 1:]
    ref_aa, alt_aa = CODON[ref_codon], CODON[alt_codon]
    if codon_i == 0 and ref_aa == "M" and alt_aa != "M":
        effect = "start_lost"
    elif ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == "*":
        effect = "stop_gained"
    elif ref_aa == "*":
        effect = "stop_lost"
    else:
        effect = "missense"
    note = ""
    if effect == "missense":
        note = "conservative" if CLASS[ref_aa] == CLASS[alt_aa] else f"{CLASS[ref_aa]} -> {CLASS[alt_aa]}"
    return {"c_change": f"c.{pos + 1}{ref}>{alt}", "codon_position": offset + 1,
            "codon_change": f"{ref_codon}>{alt_codon}", "effect": effect,
            "p_change": f"p.{THREE[ref_aa]}{codon_i + 1}{'=' if ref_aa == alt_aa else THREE[alt_aa]}",
            "type": "transition" if (ref, alt) in TRANSITIONS else "transversion", "note": note}


def annotate_indel(cds, pos, ref, alt):
    """Indel starting at 0-based pos, e.g. ref='A', alt='ATT' (insertion) or ref='ACG', alt='A' (deletion)."""
    delta = len(alt) - len(ref)
    codon_i = pos // 3
    effect = "frameshift" if delta % 3 else ("inframe_insertion" if delta > 0 else "inframe_deletion")
    new = cds[:pos] + alt + cds[pos + len(ref):]
    protein = "".join(CODON.get(new[i:i + 3], "X") for i in range(0, len(new) - 2, 3))
    stop = protein.find("*", codon_i)
    note = f"new stop {stop - codon_i} codons downstream" if effect == "frameshift" and stop >= 0 else ""
    # VCF-style anchored indels: the first base is shared, the rest is inserted or deleted.
    # (Simplified HGVS: no 3' shifting through repeats.)
    if delta > 0:
        c_change = f"c.{pos + 1}_{pos + 2}ins{alt[1:]}"
    elif len(ref) == 2:
        c_change = f"c.{pos + 2}del"
    else:
        c_change = f"c.{pos + 2}_{pos + len(ref)}del"
    return {"c_change": c_change, "codon_position": pos % 3 + 1, "codon_change": "",
            "effect": effect, "p_change": f"p.{THREE[CODON[cds[codon_i * 3:codon_i * 3 + 3]]]}{codon_i + 1}"
            + ("fs" if effect == "frameshift" else ""), "type": "indel", "note": note}


def annotate(cds, variants):
    out = []
    for pos, ref, alt in variants:
        out.append(annotate_snv(cds, pos, alt) if len(ref) == len(alt) == 1 else annotate_indel(cds, pos, ref, alt))
    return out
