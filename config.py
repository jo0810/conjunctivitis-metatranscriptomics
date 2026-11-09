"""Input and output locations shared by all Python scripts.

Run scripts from the repository root, e.g. `python 02_resistome/01_arg_quantification.py`.
  CMT_DATA_DIR    directory holding the input tables (default: data/)
  CMT_OUTPUT_DIR  directory for results (default: output/)
"""
import os

DATA_DIR = os.environ.get('CMT_DATA_DIR', 'data')
OUTPUT_DIR = os.environ.get('CMT_OUTPUT_DIR', 'output')


def in_file(*parts):
    path = os.path.join(DATA_DIR, *parts)
    if not os.path.exists(path):
        raise FileNotFoundError(f'Input file not found: {path} (see docs/input_files.md)')
    return path


def out_file(module, name):
    os.makedirs(os.path.join(OUTPUT_DIR, module), exist_ok=True)
    return os.path.join(OUTPUT_DIR, module, name)
