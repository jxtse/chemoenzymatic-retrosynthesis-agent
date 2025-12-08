#!/usr/bin/env python
"""
READRetro Helper Script
Provides convenient Python interface for READRetro biosynthesis planning
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Optional, Dict

# Add RDKit for SMILES validation
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. SMILES validation disabled.")


READRETRO_PATH = "/home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro"
CONDA_ACTIVATE = "source /opt/anaconda3/etc/profile.d/conda.sh && conda activate readretro"


def validate_smiles(smiles: str) -> Optional[str]:
    """
    Validate and canonicalize SMILES string.

    Args:
        smiles: SMILES string to validate

    Returns:
        Canonical SMILES if valid, None otherwise
    """
    if not RDKIT_AVAILABLE:
        return smiles

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return Chem.MolToSmiles(mol)


def run_single_molecule(smiles: str, gpu_id: int = 0, verbose: bool = True) -> str:
    """
    Run READRetro for a single molecule.

    Args:
        smiles: SMILES string of target molecule
        gpu_id: GPU device ID (default: 0)
        verbose: Print output to console (default: True)

    Returns:
        Output from READRetro
    """
    # Validate SMILES
    canonical = validate_smiles(smiles)
    if canonical is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Prepare command
    cmd = f"cd {READRETRO_PATH} && {CONDA_ACTIVATE} && CUDA_VISIBLE_DEVICES={gpu_id} python run.py '{canonical}'"

    if verbose:
        print(f"Running READRetro for: {canonical}")
        print(f"Using GPU: {gpu_id}")

    # Execute
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        executable='/bin/bash'
    )

    if result.returncode != 0:
        raise RuntimeError(f"READRetro failed: {result.stderr}")

    if verbose:
        print(result.stdout)

    return result.stdout


def run_batch(gpu_id: int = 0, num_threads: int = 4, verbose: bool = True) -> str:
    """
    Run READRetro in batch mode.

    Args:
        gpu_id: GPU device ID (default: 0)
        num_threads: Number of parallel threads (default: 4)
        verbose: Print output to console (default: True)

    Returns:
        Output from READRetro
    """
    cmd = f"cd {READRETRO_PATH} && {CONDA_ACTIVATE} && CUDA_VISIBLE_DEVICES={gpu_id} python run_mp.py"

    if verbose:
        print(f"Running READRetro in batch mode")
        print(f"Using GPU: {gpu_id}, Threads: {num_threads}")

    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        executable='/bin/bash'
    )

    if result.returncode != 0:
        raise RuntimeError(f"READRetro failed: {result.stderr}")

    if verbose:
        print(result.stdout)

    return result.stdout


def run_evaluation(model: str = "ensemble", gpu_id: int = 0, verbose: bool = True) -> str:
    """
    Run single-step evaluation.

    Args:
        model: Model to use ("ensemble", "retroformer", "g2s")
        gpu_id: GPU device ID (default: 0)
        verbose: Print output to console (default: True)

    Returns:
        Evaluation results
    """
    model_arg = f"-m {model}" if model != "ensemble" else ""
    cmd = f"cd {READRETRO_PATH} && {CONDA_ACTIVATE} && CUDA_VISIBLE_DEVICES={gpu_id} python eval_single.py {model_arg}"

    if verbose:
        print(f"Running evaluation with {model} model")
        print(f"Using GPU: {gpu_id}")

    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        executable='/bin/bash'
    )

    if result.returncode != 0:
        raise RuntimeError(f"Evaluation failed: {result.stderr}")

    if verbose:
        print(result.stdout)

    return result.stdout


def batch_predict(smiles_list: List[str], gpu_id: int = 0) -> Dict[str, str]:
    """
    Predict biosynthesis routes for multiple molecules.

    Args:
        smiles_list: List of SMILES strings
        gpu_id: GPU device ID (default: 0)

    Returns:
        Dictionary mapping SMILES to results
    """
    results = {}

    for i, smiles in enumerate(smiles_list):
        print(f"\nProcessing molecule {i+1}/{len(smiles_list)}: {smiles}")
        try:
            output = run_single_molecule(smiles, gpu_id=gpu_id, verbose=False)
            results[smiles] = output
            print("✓ Success")
        except Exception as e:
            print(f"✗ Failed: {e}")
            results[smiles] = f"Error: {e}"

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="READRetro Helper Script")
    parser.add_argument("--smiles", "-s", type=str, help="SMILES string")
    parser.add_argument("--file", "-f", type=str, help="File with SMILES (one per line)")
    parser.add_argument("--mode", "-m", choices=["single", "batch", "eval"], default="single")
    parser.add_argument("--gpu", "-g", type=int, default=0, help="GPU ID")
    parser.add_argument("--model", choices=["ensemble", "retroformer", "g2s"], default="ensemble")

    args = parser.parse_args()

    if args.mode == "single":
        if not args.smiles:
            print("Error: --smiles required for single mode")
            sys.exit(1)
        run_single_molecule(args.smiles, gpu_id=args.gpu)

    elif args.mode == "batch":
        if args.file:
            with open(args.file) as f:
                smiles_list = [line.strip() for line in f if line.strip()]
            results = batch_predict(smiles_list, gpu_id=args.gpu)
            print(f"\n\nProcessed {len(results)} molecules")
        else:
            run_batch(gpu_id=args.gpu)

    elif args.mode == "eval":
        run_evaluation(model=args.model, gpu_id=args.gpu)
