
import re
from collections import Counter

FIXED = {
    46: "H", 48: "H", 57: "C", 63: "H",
    71: "H", 80: "H", 83: "D", 120: "H",
    124: "D", 143: "R", 146: "C"
}


def parse_fasta(text):
    records = []
    header = None
    chunks = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks)))

    return records


def get_score(header, key):
    match = re.search(
        rf"(?:^|,\s*){re.escape(key)}=([0-9.eE+-]+)",
        header
    )
    return float(match.group(1)) if match else None


def analyze_fasta(text):
    records = parse_fasta(text)

    if not records:
        raise ValueError("FASTAに配列がありません")

    reference = records[0][1].split("/")[0]

    if len(reference) != 153:
        raise ValueError("参照配列が153残基ではありません")

    designs = [
        (header, sequence)
        for header, sequence in records[1:]
        if header.startswith("T=")
    ]

    counts = Counter(
        sequence.split("/")[0]
        for _, sequence in designs
    )

    results = []

    for number, (header, sequence) in enumerate(designs, 1):
        chains = sequence.split("/")
        valid = (
            len(chains) == 2
            and all(len(chain) == 153 for chain in chains)
            and chains[0] == chains[1]
        )

        if valid:
            valid = all(
                chain[position - 1] == amino_acid
                for chain in chains
                for position, amino_acid in FIXED.items()
            )

        chain_a = chains[0]

        mutations = [
            f"{original}{position}{changed}"
            for position, (original, changed)
            in enumerate(zip(reference, chain_a), 1)
            if original != changed
        ] if len(chain_a) == len(reference) else []

        results.append({
            "number": number,
            "sequence": chain_a,
            "valid": valid,
            "mutation_count": len(mutations),
            "mutations": ", ".join(mutations),
            "score": get_score(header, "score"),
            "global_score": get_score(header, "global_score"),
            "recovery": get_score(header, "seq_recovery"),
            "duplicate_count": counts[chain_a]
        })

    return {
        "reference": reference,
        "results": results,
        "total": len(results),
        "valid": sum(r["valid"] for r in results),
        "unique": len(counts)
    }
