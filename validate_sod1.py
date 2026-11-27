import subprocess

FASTA = (
    "/home/kyokugen/protein-work/ProteinMPNN/"
    "sod1_project/results/gui_20261120_151930/seqs/2C9V.fa"
)

FIXED = {
    46: "H",
    48: "H",
    57: "C",
    63: "H",
    71: "H",
    80: "H",
    83: "D",
    120: "H",
    124: "D",
    143: "R",
    146: "C",
}

result = subprocess.run(
    ["wsl.exe", "-d", "Ubuntu", "--", "cat", FASTA],
    capture_output=True,
    text=True,
    encoding="utf-8",
    check=True
)

records = []
header = None
parts = []

for line in result.stdout.splitlines():
    line = line.strip()

    if line.startswith(">"):
        if header is not None:
            records.append((header, "".join(parts)))
        header = line
        parts = []
    elif line:
        parts.append(line)

if header is not None:
    records.append((header, "".join(parts)))

designed = [
    sequence
    for header, sequence in records
    if header.startswith(">T=")
]

passed = 0

for i, sequence in enumerate(designed, 1):
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

    if valid:
        passed += 1

    print(f"Sequence {i}: {'PASS' if valid else 'FAIL'}")

print(f"Fixed residue validation: {passed}/{len(designed)} PASS")