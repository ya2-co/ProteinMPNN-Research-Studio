
import subprocess
from sod1_analysis import analyze_fasta

path = (
    "/home/kyokugen/protein-work/ProteinMPNN/"
    "sod1_project/results/gui_20261120_153239/seqs/2C9V.fa"
)

fasta = subprocess.check_output(
    ["wsl.exe", "-d", "Ubuntu", "--", "cat", path],
    text=True,
    encoding="utf-8"
)

data = analyze_fasta(fasta)

print("生成配列数:", data["total"])
print("制約検証:", data["valid"], "/", data["total"])
print("ユニーク配列数:", data["unique"])

for r in data["results"][:3]:
    print(
        "配列:", r["number"],
        "変異数:", r["mutation_count"],
        "Score:", r["score"],
        "Recovery:", r["recovery"],
        "重複数:", r["duplicate_count"]
    )
