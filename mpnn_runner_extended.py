import subprocess
from datetime import datetime
from pathlib import Path

BASE = '/home/kyokugen/protein-work/ProteinMPNN'
PROJECT = f'{BASE}/sod1_project'
PYTHON = '/home/kyokugen/protein-work/.venv/bin/python'


def run_mpnn(num_sequences=10, temperature=0.1, seed=42, *,
             model_name='v_48_020', batch_size=1, backbone_noise=0.0,
             use_soluble_model=False, ca_only=False):
    """Run official ProteinMPNN with the existing SOD1 constraints.

    This keeps the existing GUI's three positional arguments compatible.
    """
    if isinstance(num_sequences, bool) or not isinstance(num_sequences, int) or num_sequences < 1:
        raise ValueError('num_sequences must be a positive integer')
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError('batch_size must be a positive integer')
    if num_sequences % batch_size:
        raise ValueError('num_sequences must be divisible by batch_size')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    if not (0 < float(temperature) <= 1):
        raise ValueError('temperature must be between 0 and 1')
    if float(backbone_noise) < 0:
        raise ValueError('backbone_noise must be nonnegative')
    if ca_only and use_soluble_model:
        raise ValueError('ca_only and use_soluble_model cannot be combined')
    if not model_name or any(ch not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for ch in model_name):
        raise ValueError('invalid model_name')

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output = f'{PROJECT}/results/gui_{timestamp}'
    command = [
        'wsl.exe', '-d', 'Ubuntu', '--', PYTHON,
        f'{BASE}/protein_mpnn_run.py',
        '--jsonl_path', f'{PROJECT}/config/parsed.jsonl',
        '--chain_id_jsonl', f'{PROJECT}/config/chains.jsonl',
        '--fixed_positions_jsonl', f'{PROJECT}/config/fixed.jsonl',
        '--tied_positions_jsonl', f'{PROJECT}/config/tied.jsonl',
        '--out_folder', output,
        '--num_seq_per_target', str(num_sequences),
        '--sampling_temp', str(temperature),
        '--seed', str(seed),
        '--batch_size', str(batch_size),
        '--model_name', model_name,
        '--backbone_noise', str(backbone_noise),
    ]
    if use_soluble_model:
        command.extend(['--use_soluble_model', '1'])
    if ca_only:
        command.extend(['--ca_only', '1'])

    print('Starting ProteinMPNN...')
    print('Output:', output)
    result = subprocess.run(command, capture_output=True, text=True,
                            encoding='utf-8', errors='replace')
    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError('ProteinMPNN execution failed:\n' + result.stderr[-4000:])
    print('ProteinMPNN completed successfully.')
    return output


if __name__ == '__main__':
    run_mpnn()
