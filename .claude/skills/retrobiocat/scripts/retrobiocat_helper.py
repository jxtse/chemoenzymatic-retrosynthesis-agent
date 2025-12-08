#!/usr/bin/env python
"""
RetroBioCat 2.0 Helper Script
Provides convenient interface for synthesis planning
"""

import argparse
import json
from typing import List, Dict, Optional
from pathlib import Path


def single_step_retrosynthesis(target_smi: str,
                               expanders: List[str] = None,
                               verbose: bool = True) -> List:
    """
    Perform single-step retrosynthesis.

    Args:
        target_smi: Target molecule SMILES
        expanders: List of expander names (default: ['retrobiocat'])
        verbose: Print detailed output

    Returns:
        List of Reaction objects
    """
    from rbc2 import get_expanders

    if expanders is None:
        expanders = ['retrobiocat']

    if verbose:
        print(f"Target: {target_smi}")
        print(f"Using expanders: {', '.join(expanders)}")

    expanders_obj = get_expanders(expanders)

    all_reactions = []
    for expander in expanders_obj:
        reactions = expander.get_reactions(target_smi)
        all_reactions.extend(reactions)

        if verbose:
            print(f"\n{expander.__class__.__name__}: {len(reactions)} reactions")

    if verbose:
        print(f"\nTotal reactions found: {len(all_reactions)}")

        # Show top reactions
        for i, rxn in enumerate(all_reactions[:5], 1):
            print(f"\n{i}. {rxn.reaction_smiles()}")
            print(f"   Score: {rxn.score:.3f}")
            print(f"   Type: {rxn.rxn_type}")
            print(f"   Domain: {rxn.rxn_domain}")

            if rxn.precedents:
                print(f"   Precedents: {len(rxn.precedents)}")

    return all_reactions


def multi_step_retrosynthesis(target_smi: str,
                              expanders: List[str] = None,
                              max_time: int = 30,
                              max_depth: int = 10,
                              use_commercial_sme: bool = True,
                              verbose: bool = True) -> List:
    """
    Perform multi-step retrosynthesis using MCTS.

    Args:
        target_smi: Target molecule SMILES
        expanders: List of expander names
        max_time: Maximum search time (seconds)
        max_depth: Maximum pathway depth
        use_commercial_sme: Use commercial availability checker
        verbose: Print detailed output

    Returns:
        List of solved Pathway objects
    """
    from rbc2 import MCTS, get_expanders, CommercialSME

    if expanders is None:
        expanders = ['retrobiocat', 'aizynthfinder']

    if verbose:
        print(f"Target: {target_smi}")
        print(f"Expanders: {', '.join(expanders)}")
        print(f"Max time: {max_time}s, Max depth: {max_depth}")

    # Initialize expanders
    expanders_obj = get_expanders(expanders)

    # Initialize SME
    sme = CommercialSME() if use_commercial_sme else None

    # Run MCTS
    mcts = MCTS(target_smi, expanders_obj, starting_material_evaluator=sme)
    mcts.config.max_search_time = max_time
    mcts.config.max_depth = max_depth

    if verbose:
        print("\nRunning MCTS search...")

    mcts.run()

    # Get results
    all_pathways = mcts.get_all_pathways()
    solved_pathways = mcts.get_solved_pathways()

    if verbose:
        print(f"\nSearch complete!")
        print(f"Total pathways: {len(all_pathways)}")
        print(f"Solved pathways: {len(solved_pathways)}")

        # Show top pathways
        for i, pathway in enumerate(solved_pathways[:3], 1):
            print(f"\n=== Pathway {i} ===")
            print(f"Length: {pathway.pathway_length} steps")
            print(f"Reactions: {len(pathway.reactions)}")

            bio_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'biocatalysis')
            chem_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'chemistry')
            print(f"Bio steps: {bio_rxns}, Chem steps: {chem_rxns}")

            print("\nReaction sequence:")
            for j, rxn in enumerate(pathway.reactions, 1):
                print(f"  {j}. {rxn.reaction_smiles()} ({rxn.rxn_domain})")

    return solved_pathways


def check_availability(smiles_list: List[str],
                      evaluator: str = 'commercial',
                      verbose: bool = True) -> Dict[str, tuple]:
    """
    Check starting material availability.

    Args:
        smiles_list: List of SMILES to check
        evaluator: 'commercial' or 'ecoli'
        verbose: Print results

    Returns:
        Dictionary mapping SMILES to (available, info) tuples
    """
    from rbc2 import CommercialSME, EcoliSME

    if evaluator == 'commercial':
        sme = CommercialSME()
    elif evaluator == 'ecoli':
        sme = EcoliSME()
    else:
        raise ValueError(f"Unknown evaluator: {evaluator}")

    results = {}

    if verbose:
        print(f"Checking availability using {evaluator} evaluator\n")

    for smi in smiles_list:
        available, info = sme.eval(smi)
        results[smi] = (available, info)

        if verbose:
            print(f"{smi}")
            print(f"  Available: {available}")
            if info:
                print(f"  Info: {info}")
            print()

    return results


def batch_planning(smiles_file: str,
                  expanders: List[str] = None,
                  max_time: int = 30,
                  output_file: str = None,
                  verbose: bool = True) -> Dict[str, List]:
    """
    Plan routes for multiple targets from file.

    Args:
        smiles_file: File with SMILES (one per line)
        expanders: List of expander names
        max_time: Max time per target
        output_file: Output JSON file
        verbose: Print progress

    Returns:
        Dictionary mapping SMILES to pathways
    """
    from rbc2 import MCTS, get_expanders

    # Read targets
    with open(smiles_file) as f:
        targets = [line.strip() for line in f if line.strip()]

    if expanders is None:
        expanders = ['retrobiocat', 'aizynthfinder']

    expanders_obj = get_expanders(expanders)

    results = {}

    for i, target in enumerate(targets, 1):
        if verbose:
            print(f"\n[{i}/{len(targets)}] Processing {target}")

        try:
            mcts = MCTS(target, expanders_obj)
            mcts.config.max_search_time = max_time
            mcts.run()

            pathways = mcts.get_solved_pathways()
            results[target] = pathways

            if verbose:
                print(f"  Found {len(pathways)} pathways")

        except Exception as e:
            if verbose:
                print(f"  Error: {e}")
            results[target] = []

    # Save if requested
    if output_file:
        output_data = {}
        for target, pathways in results.items():
            output_data[target] = [p.save() for p in pathways]

        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)

        if verbose:
            print(f"\nResults saved to {output_file}")

    return results


def analyze_pathway(pathway, sme=None, verbose: bool = True):
    """
    Analyze a pathway in detail.

    Args:
        pathway: Pathway object
        sme: Starting material evaluator (optional)
        verbose: Print analysis
    """
    if verbose:
        print("=== Pathway Analysis ===")
        print(f"Target: {pathway.target_smi}")
        print(f"Length: {pathway.pathway_length} steps")
        print(f"Total reactions: {len(pathway.reactions)}")
        print(f"Total molecules: {len(pathway.all_smis)}")

        # Domain breakdown
        bio_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'biocatalysis')
        chem_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'chemistry')
        biosyn_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'biosynthesis')

        print(f"\nDomain breakdown:")
        print(f"  Biocatalysis: {bio_rxns}")
        print(f"  Chemistry: {chem_rxns}")
        print(f"  Biosynthesis: {biosyn_rxns}")

        # Starting materials
        end_smis = pathway.end_smis()
        print(f"\nStarting materials ({len(end_smis)}):")
        for smi in end_smis:
            print(f"  {smi}")

            if sme:
                available, info = sme.eval(smi)
                print(f"    Available: {available}")

        # Reaction details
        print("\nReactions:")
        for i, rxn in enumerate(pathway.reactions, 1):
            print(f"\n  {i}. {rxn.name}")
            print(f"     {rxn.reaction_smiles()}")
            print(f"     Type: {rxn.rxn_type}, Domain: {rxn.rxn_domain}")
            print(f"     Score: {rxn.score:.3f}")

            if rxn.precedents:
                print(f"     Precedents: {len(rxn.precedents)}")
                top_prec = rxn.precedents[0]
                print(f"       Top: {top_prec.name} (sim: {top_prec.similarity:.3f})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="RetroBioCat 2.0 Helper Script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single-step retrosynthesis
  python retrobiocat_helper.py single -s "CCCC=O" -e retrobiocat enzymemap

  # Multi-step planning
  python retrobiocat_helper.py multi -s "CC(O)C(=O)O" -t 60 -d 12

  # Check availability
  python retrobiocat_helper.py check -s "CCO" "c1ccccc1O" -e commercial

  # Batch processing
  python retrobiocat_helper.py batch -f targets.txt -o results.json
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Single-step command
    single_parser = subparsers.add_parser('single', help='Single-step retrosynthesis')
    single_parser.add_argument('-s', '--smiles', required=True, help='Target SMILES')
    single_parser.add_argument('-e', '--expanders', nargs='+',
                              default=['retrobiocat'],
                              help='Expanders to use')

    # Multi-step command
    multi_parser = subparsers.add_parser('multi', help='Multi-step retrosynthesis')
    multi_parser.add_argument('-s', '--smiles', required=True, help='Target SMILES')
    multi_parser.add_argument('-e', '--expanders', nargs='+',
                             default=['retrobiocat', 'aizynthfinder'],
                             help='Expanders to use')
    multi_parser.add_argument('-t', '--time', type=int, default=30,
                             help='Max search time (seconds)')
    multi_parser.add_argument('-d', '--depth', type=int, default=10,
                             help='Max pathway depth')
    multi_parser.add_argument('--no-commercial', action='store_true',
                             help='Disable commercial availability check')

    # Check command
    check_parser = subparsers.add_parser('check', help='Check availability')
    check_parser.add_argument('-s', '--smiles', nargs='+', required=True,
                             help='SMILES to check')
    check_parser.add_argument('-e', '--evaluator', choices=['commercial', 'ecoli'],
                             default='commercial', help='Evaluator to use')

    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Batch processing')
    batch_parser.add_argument('-f', '--file', required=True,
                             help='File with SMILES (one per line)')
    batch_parser.add_argument('-e', '--expanders', nargs='+',
                             default=['retrobiocat', 'aizynthfinder'],
                             help='Expanders to use')
    batch_parser.add_argument('-t', '--time', type=int, default=30,
                             help='Max time per target (seconds)')
    batch_parser.add_argument('-o', '--output', help='Output JSON file')

    args = parser.parse_args()

    if args.command == 'single':
        single_step_retrosynthesis(args.smiles, args.expanders)

    elif args.command == 'multi':
        multi_step_retrosynthesis(
            args.smiles,
            args.expanders,
            args.time,
            args.depth,
            not args.no_commercial
        )

    elif args.command == 'check':
        check_availability(args.smiles, args.evaluator)

    elif args.command == 'batch':
        batch_planning(args.file, args.expanders, args.time, args.output)

    else:
        parser.print_help()
