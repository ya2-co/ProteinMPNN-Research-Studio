import subprocess
from datetime import datetime

BASE = "/home/kyokugen/protein-work/ProteinMPNN"
PROJECT = f"{BASE}/sod1_project"
PYTHON = "/home/kyokugen/protein-work/.venv/bin/python"


def run_mpnn(num_sequences=10, temperature=0.1, seed=42):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = f"{PROJECT}/results/gui_{timestamp}"

    command = [
        "wsl.exe", "-d", "Ubuntu", "--",
        PYTHON,
        f"{BASE}/protein_mpnn_run.py",
        "--jsonl_path", f"{PROJECT}/config/parsed.jsonl",
        "--chain_id_jsonl", f"{PROJECT}/config/chains.jsonl",
        "--fixed_positions_jsonl", f"{PROJECT}/config/fixed.jsonl",
        "--tied_positions_jsonl", f"{PROJECT}/config/tied.jsonl",
        "--out_folder", output,
        "--num_seq_per_target", str(num_sequences),
        "--sampling_temp", str(temperature),
        "--seed", str(seed),
        "--batch_size", "1",
        "--model_name", "v_48_020",
    ]

    print("Starting ProteinMPNN...")
    print("Output:", output)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace"
    )

    print(result.stdout)

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("ProteinMPNN execution failed")

    print("ProteinMPNN completed successfully.")
    return output


if __name__ == "__main__":
    run_mpnn()