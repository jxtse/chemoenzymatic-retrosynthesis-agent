---
name: retrobiocat
description: "Modular biocatalysis retrosynthesis framework. Combine multiple approaches (RetroBioCat, EnzymeMap, BKMS, RetroRules, AIZynthfinder) for hybrid chemistry+biocatalysis synthesis planning with MCTS, literature precedents, commercial availability checks."
---

# RetroBioCat 2.0: Computer-Aided Synthesis Planning

## Overview

RetroBioCat 2.0 (rbc2) is a modular Python package for computer-aided synthesis planning in biocatalysis. It provides a flexible framework for combining multiple retrosynthesis approaches to enable hybrid synthesis planning that integrates both biocatalysis and traditional chemistry.

**Core Capabilities:**
- Multiple retrosynthesis methods (expanders) in one framework
- Hybrid chemistry + biocatalysis route planning
- Monte Carlo Tree Search (MCTS) for multi-step synthesis
- Literature precedent searching (BRENDA, BKMS, custom databases)
- Commercial availability and E. coli metabolism evaluation
- Modular, extensible architecture for custom expanders

**Key Distinction:** RetroBioCat 2.0 is a **Python package only** (no web interface) designed for maximum flexibility and integration with custom workflows.

## When to Use This Skill

This skill should be used when:

- "Plan biocatalytic synthesis route"
- "Combine enzymatic and chemical steps"
- "Compare multiple retrosynthesis approaches"
- "Find enzymatic reactions for this transformation"
- "Hybrid synthesis planning with enzymes and chemistry"
- "Check commercial availability of precursors"
- "Multi-step biocatalysis route design"
- Need to integrate BRENDA, BKMS, or RetroRules data
- Want literature precedents for enzymatic reactions

## Installation and Environment Setup

### Installation with uv

RetroBioCat 2.0 is installed via pip from GitHub. We use uv for package management:

```bash
# Install into current uv project
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

**Requirements:**
- Python >= 3.9 (compatible with Python 3.12)
- RDKit (installed automatically)
- Additional dependencies installed automatically on first use

### ARM Mac Installation

For Apple Silicon Macs, install HDF5 first:

```bash
brew install hdf5 c-blosc
export HDF5_DIR=/opt/homebrew/opt/hdf5
export BLOSC_DIR=/opt/homebrew/opt/c-blosc
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### Verify Installation

```bash
# Verify using uv
uv run python -c "import rbc2; from rbc2 import RetroBioCatExpander; print('RetroBioCat 2.0 installed successfully!')"
```

Or in Python:

```python
import rbc2
from rbc2 import RetroBioCatExpander

print("RetroBioCat 2.0 installed successfully!")

# Test basic functionality
expander = RetroBioCatExpander()
reactions = expander.get_reactions('CCCC=O')
print(f"Found {len(reactions)} reactions")
```

## Available Expanders

### Biocatalysis Expanders

#### 1. RetroBioCat
Manually curated enzyme reaction rules with literature precedent search.

```python
from rbc2 import RetroBioCatExpander

expander = RetroBioCatExpander(
    include_experimental=False,
    include_two_step=True,
    include_requires_absence_of_water=False,
    score_similarity_before_option_creation=True,
    search_literature_precedent=True,
    only_active_literature_precedent=True,
    similarity_cutoff=0.55
)

reactions = expander.get_reactions('CCCC=O')
```

**Use for:** Well-characterized enzymatic transformations, biocatalysis cascades

**Reference:** Finnigan, W., et al. Nature Catalysis 4, 98–104 (2021)

#### 2. EnzymeMap
BRENDA-based templates with automated extraction and relevance scoring.

```python
from rbc2 import EnzymeMapExpander

expander = EnzymeMapExpander(
    cutoff_cumulative=0.995,
    cutoff_number=50,
    enable_precedent_search=True,
    similarity_cutoff=0.1
)
```

**Use for:** Broad enzyme coverage, novel enzymatic transformations

**Reference:** Heid, E., et al. Chem. Sci. 14, 14229–14242 (2023)

#### 3. BKMS
BKMS database templates with similarity-based precedent search.

```python
from rbc2 import BKMSExpander

expander = BKMSExpander(
    cutoff_cumulative=0.995,
    cutoff_number=50,
    allow_multi_product_templates=False,
    enable_precedent_search=True,
    similarity_cutoff=0.1
)
```

**Use for:** Metabolic pathways, biosynthesis routes

**Reference:** Levin, I., et al. Nat Commun 13, 7747 (2022)

#### 4. RetroRules
RetroPath templates with biological scoring.

```python
from rbc2 import RetroRulesExpander

expander = RetroRulesExpander(
    rank_by='combined_score',  # options: similarity, score, combined_score
    diameter=6,
    similarity_threshold=0.2,
    score_threshold=0.2,
    combined_score_threshold=0.2,
    max_reactions=100
)
```

**Use for:** Metabolic engineering, pathway design

**Reference:** Koch, M., et al. ACS Synth. Biol. 9, 157–168 (2020)

### Chemistry Expanders

#### 5. AIZynthfinder
AI-based retrosynthesis for chemical transformations.

```python
from rbc2 import AIZynthfinderExpander

expander = AIZynthfinderExpander(
    cutoff_cumulative=0.995,
    cutoff_number=50
)
```

**Use for:** Traditional organic chemistry, hybrid routes

**Reference:** Genheden, S., et al. J Cheminform 12, 70 (2020)

#### 6. RingBreaker
Neural network for ring system synthesis.

```python
from rbc2 import RingBreakerPolicyExpander

expander = RingBreakerPolicyExpander(
    cutoff_cumulative=0.995,
    cutoff_number=10
)
```

**Use for:** Complex ring formations, cyclization reactions

**Reference:** Thakkar, A., et al. J. Med. Chem. 63, 8791–8808 (2020)

#### 7. AskCos
Robotic platform-informed synthesis planning.

```python
from rbc2 import AskcosPolicyExpander

expander = AskcosPolicyExpander(
    cutoff_cumulative=0.995,
    cutoff_number=50
)
```

**Use for:** Automated synthesis, flow chemistry

**Reference:** Coley, C. W., et al. Science 365, eaax1566 (2019)

## Core Workflows

### Workflow 1: Single-Step Retrosynthesis

**Use Case:** Explore one-step synthetic options for a target molecule

**Example - Using RetroBioCat:**

```python
from rbc2 import RetroBioCatExpander

# Initialize expander
expander = RetroBioCatExpander()

# Get reactions for target molecule (butanal)
reactions = expander.get_reactions('CCCC=O')

# Explore reactions
for rxn in reactions:
    print(f"Reaction: {rxn.reaction_smiles()}")
    print(f"Score: {rxn.score}")
    print(f"Name: {rxn.name}")
    print(f"Type: {rxn.rxn_type}")

    # Check literature precedents
    for precedent in rxn.precedents:
        print(f"  Precedent: {precedent.name}")
        print(f"  Similarity: {precedent.similarity}")
        print(f"  Data: {precedent.data}")

    # Complexity change
    complexity_change = rxn.get_complexity_change()
    print(f"Complexity change: {complexity_change}\n")
```

**Example - Comparing Multiple Expanders:**

```python
from rbc2 import RetroBioCatExpander, EnzymeMapExpander, AIZynthfinderExpander

target = 'CC(O)C(=O)O'  # Lactic acid

expanders = {
    'RetroBioCat': RetroBioCatExpander(),
    'EnzymeMap': EnzymeMapExpander(),
    'AIZynthfinder': AIZynthfinderExpander()
}

for name, expander in expanders.items():
    reactions = expander.get_reactions(target)
    print(f"{name}: Found {len(reactions)} reactions")

    # Show top reaction
    if reactions:
        top_rxn = reactions[0]
        print(f"  Top: {top_rxn.reaction_smiles()}")
        print(f"  Score: {top_rxn.score}\n")
```

### Workflow 2: Multi-Step Retrosynthesis with MCTS

**Use Case:** Plan complete synthesis routes using Monte Carlo Tree Search

**Basic Example:**

```python
from rbc2 import MCTS, get_expanders

# Target molecule
target_smi = 'CC(O)C(=O)O'  # Lactic acid

# Get expanders (convenience function)
expanders = get_expanders(['retrobiocat', 'aizynthfinder'])

# Initialize MCTS
mcts = MCTS(target_smi, expanders)

# Configure search
mcts.config.max_search_time = 30  # seconds
mcts.config.max_iterations = 1000
mcts.config.c_exploration = 1.414  # UCB exploration constant

# Run search
mcts.run()

# Get results
all_pathways = mcts.get_all_pathways()
solved_pathways = mcts.get_solved_pathways()

print(f"Found {len(all_pathways)} total pathways")
print(f"Found {len(solved_pathways)} solved pathways")

# Analyze best pathway
if solved_pathways:
    best_pathway = solved_pathways[0]
    print(f"\nBest pathway length: {best_pathway.pathway_length}")
    print(f"Reactions in pathway: {len(best_pathway.reactions)}")

    for i, rxn in enumerate(best_pathway.reactions):
        print(f"Step {i+1}: {rxn.reaction_smiles()}")
```

**Advanced Example with Custom Configuration:**

```python
from rbc2 import MCTS, RetroBioCatExpander, EnzymeMapExpander, AIZynthfinderExpander
from rbc2 import CommercialSME

# Target molecule
target_smi = 'c1ccc2c(c1)c(c(=O)oc2=O)O'  # 4-Hydroxycoumarin

# Initialize expanders
expanders = [
    RetroBioCatExpander(similarity_cutoff=0.6),
    EnzymeMapExpander(cutoff_number=30),
    AIZynthfinderExpander()
]

# Initialize starting material evaluator
sme = CommercialSME()

# Initialize MCTS
mcts = MCTS(target_smi, expanders, starting_material_evaluator=sme)

# Configure
mcts.config.max_search_time = 60
mcts.config.max_depth = 10
mcts.config.c_exploration = 1.5

# Run
mcts.run()

# Get solved pathways
solved_pathways = mcts.get_solved_pathways()

# Sort by pathway length
solved_pathways.sort(key=lambda p: p.pathway_length)

# Display results
for i, pathway in enumerate(solved_pathways[:5]):  # Top 5
    print(f"\n=== Pathway {i+1} ===")
    print(f"Length: {pathway.pathway_length} steps")
    print(f"Reactions: {len(pathway.reactions)}")

    # Show reaction sequence
    for j, rxn in enumerate(pathway.reactions):
        print(f"  {j+1}. {rxn.name}")
        print(f"     {rxn.reaction_smiles()}")
        print(f"     Domain: {rxn.rxn_domain}")
```

### Workflow 3: Starting Material Evaluation

**Use Case:** Check commercial availability or metabolic presence

**Commercial Availability:**

```python
from rbc2 import CommercialSME

sme = CommercialSME()

molecules = [
    'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
    'CC(C)CC(C)O',             # Isopentanol
    'c1ccccc1O'                # Phenol
]

for smi in molecules:
    available, info = sme.eval(smi)
    print(f"{smi}")
    print(f"  Available: {available}")
    print(f"  Info: {info}\n")
```

**E. coli Metabolism:**

```python
from rbc2 import EcoliSME

sme = EcoliSME()

metabolites = [
    'C(C(=O)O)N',           # Glycine
    'CC(C(=O)O)N',          # Alanine
    'C1=CC=C(C=C1)C(=O)O'  # Benzoic acid
]

for smi in metabolites:
    in_ecoli, info = sme.eval(smi)
    print(f"{smi}")
    print(f"  In E. coli: {in_ecoli}")
    print(f"  Info: {info}\n")
```

**Custom Starting Material Evaluator:**

```python
from rbc2 import StartingMaterialEvaluatorInterface
from typing import Tuple
from rdkit import Chem

class CustomSME(StartingMaterialEvaluatorInterface):
    """Custom evaluator for lab-specific building blocks"""

    def __init__(self, custom_molecules: list):
        super().__init__()
        self.custom_molecules = set(custom_molecules)

    def eval(self, smi: str) -> Tuple[bool, dict]:
        # Canonicalize
        mol = Chem.MolFromSmiles(smi)
        canonical = Chem.MolToSmiles(mol)

        available = canonical in self.custom_molecules
        info = {'source': 'custom_lab_inventory'}

        return available, info

    def is_mol_chiral(self, smi: str) -> bool:
        mol = Chem.MolFromSmiles(smi)
        return Chem.FindMolChiralCenters(mol, includeUnassigned=True)

# Usage
lab_inventory = ['CC(=O)O', 'CCO', 'c1ccccc1']
custom_sme = CustomSME(lab_inventory)

# Use in MCTS
mcts = MCTS(target_smi, expanders, starting_material_evaluator=custom_sme)
```

### Workflow 4: Pathway Analysis and Export

**Use Case:** Analyze and export planned pathways

```python
from rbc2 import MCTS, get_expanders, CommercialSME, load_pathway
import json

# Run MCTS
target_smi = 'CC(O)C(=O)O'
expanders = get_expanders(['retrobiocat', 'enzymemap'])
sme = CommercialSME()

mcts = MCTS(target_smi, expanders, starting_material_evaluator=sme)
mcts.config.max_search_time = 30
mcts.run()

solved_pathways = mcts.get_solved_pathways()

if solved_pathways:
    pathway = solved_pathways[0]

    # Basic pathway info
    print(f"Target: {pathway.target_smi}")
    print(f"Length: {pathway.pathway_length}")
    print(f"Total reactions: {len(pathway.reactions)}")

    # Get end molecules (starting materials)
    end_smiles = pathway.end_smis()
    print(f"\nStarting materials: {end_smiles}")

    # Check if all are commercially available
    for smi in end_smiles:
        available, info = sme.eval(smi)
        print(f"  {smi}: {available}")

    # Get all molecules in pathway
    print(f"\nAll molecules: {len(pathway.all_smis)}")

    # Get tree structure
    tree = pathway.tree
    print(f"\nTree structure: {json.dumps(tree, indent=2)}")

    # Export pathway
    pathway_dict = pathway.save()

    # Save to file
    with open('pathway_output.json', 'w') as f:
        json.dump(pathway_dict, f, indent=2)

    # Load pathway later
    loaded_pathway = load_pathway(pathway_dict)

    # Get PARoutes format
    pa_routes = pathway.get_pa_route(sme)
    print(f"\nPARoutes format: {json.dumps(pa_routes, indent=2)}")
```

## Advanced Usage

### Combining Multiple Expanders Strategically

```python
from rbc2 import (RetroBioCatExpander, EnzymeMapExpander,
                  AIZynthfinderExpander, MCTS)

# Strategy: Try biocatalysis first, fall back to chemistry
bio_expanders = [
    RetroBioCatExpander(similarity_cutoff=0.7),  # Strict
    EnzymeMapExpander(cutoff_number=20)
]

chem_expanders = [
    AIZynthfinderExpander()
]

# First try bio-only
mcts_bio = MCTS(target_smi, bio_expanders)
mcts_bio.config.max_search_time = 15
mcts_bio.run()

bio_pathways = mcts_bio.get_solved_pathways()

if not bio_pathways:
    print("No bio-only routes found, trying hybrid...")
    # Try hybrid
    mcts_hybrid = MCTS(target_smi, bio_expanders + chem_expanders)
    mcts_hybrid.config.max_search_time = 30
    mcts_hybrid.run()

    hybrid_pathways = mcts_hybrid.get_solved_pathways()
    print(f"Found {len(hybrid_pathways)} hybrid pathways")
else:
    print(f"Found {len(bio_pathways)} bio-only pathways")
```

### Custom Expander Configuration

```python
from rbc2.configs import Expansion_Config
from rbc2 import RetroBioCatExpander

# Create custom config
config = Expansion_Config()
config.max_reactions = 20           # Max reactions per step
config.min_reaction_score = 0.3     # Score threshold
config.return_full_precedents = True
config.max_time = 60                # Max time for expansion

# Apply to expander
expander = RetroBioCatExpander(config=config)

# Or modify after initialization
expander.config.max_reactions = 15
```

### Batch Processing Multiple Targets

```python
from rbc2 import MCTS, get_expanders, CommercialSME
import pandas as pd

# List of targets
targets = [
    ('Lactic acid', 'CC(O)C(=O)O'),
    ('Phenylethanol', 'c1ccc(cc1)CCO'),
    ('Vanillin', 'COc1cc(ccc1O)C=O')
]

# Setup
expanders = get_expanders(['retrobiocat', 'enzymemap', 'aizynthfinder'])
sme = CommercialSME()

results = []

for name, smiles in targets:
    print(f"\nProcessing {name}...")

    mcts = MCTS(smiles, expanders, starting_material_evaluator=sme)
    mcts.config.max_search_time = 20
    mcts.run()

    solved = mcts.get_solved_pathways()

    results.append({
        'name': name,
        'smiles': smiles,
        'pathways_found': len(solved),
        'min_length': min([p.pathway_length for p in solved]) if solved else None,
        'avg_length': sum([p.pathway_length for p in solved]) / len(solved) if solved else None
    })

# Create summary DataFrame
df = pd.DataFrame(results)
print("\n=== Summary ===")
print(df)
```

## Common Use Cases and Examples

### Use Case 1: Natural Product Synthesis

```python
from rbc2 import MCTS, get_expanders

# Target: Vanillin (natural product)
target = 'COc1cc(ccc1O)C=O'

# Use biocatalysis-focused expanders
expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms'])

mcts = MCTS(target, expanders)
mcts.config.max_search_time = 45
mcts.run()

pathways = mcts.get_solved_pathways()

# Filter for bio-only routes
bio_only = [p for p in pathways
            if all(rxn.rxn_domain == 'biocatalysis' for rxn in p.reactions)]

print(f"Found {len(bio_only)} biocatalysis-only routes")
```

### Use Case 2: Pharmaceutical Intermediate

```python
from rbc2 import MCTS, get_expanders

# Target: Ibuprofen intermediate
target = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'

# Hybrid approach
expanders = get_expanders(['retrobiocat', 'aizynthfinder', 'askcos'])

mcts = MCTS(target, expanders)
mcts.config.max_search_time = 60
mcts.config.max_depth = 8
mcts.run()

pathways = mcts.get_solved_pathways()

# Analyze domain distribution
for i, pathway in enumerate(pathways[:3]):
    bio_steps = sum(1 for rxn in pathway.reactions if rxn.rxn_domain == 'biocatalysis')
    chem_steps = sum(1 for rxn in pathway.reactions if rxn.rxn_domain == 'chemistry')

    print(f"Pathway {i+1}: {bio_steps} bio, {chem_steps} chem steps")
```

### Use Case 3: Metabolic Engineering

```python
from rbc2 import MCTS, get_expanders, EcoliSME

# Target: Non-natural amino acid
target = 'CC(C)C(N)C(=O)O'  # Valine

# Use E. coli as host
sme = EcoliSME()
expanders = get_expanders(['retrorules', 'bkms', 'enzymemap'])

mcts = MCTS(target, expanders, starting_material_evaluator=sme)
mcts.config.max_search_time = 30
mcts.run()

pathways = mcts.get_solved_pathways()

# Filter pathways that start from E. coli metabolites
ecoli_pathways = []
for pathway in pathways:
    end_smis = pathway.end_smis()
    all_in_ecoli = all(sme.eval(smi)[0] for smi in end_smis)

    if all_in_ecoli:
        ecoli_pathways.append(pathway)

print(f"Found {len(ecoli_pathways)} pathways from E. coli metabolites")
```

### Use Case 4: Literature Precedent Analysis

```python
from rbc2 import RetroBioCatExpander, EnzymeMapExpander

target = 'c1ccc(cc1)CCO'  # Phenylethanol

# Get reactions with precedents
rbc_expander = RetroBioCatExpander(search_literature_precedent=True)
em_expander = EnzymeMapExpander(enable_precedent_search=True)

rbc_reactions = rbc_expander.get_reactions(target)
em_reactions = em_expander.get_reactions(target)

# Analyze precedents
print("=== RetroBioCat Precedents ===")
for rxn in rbc_reactions[:3]:
    print(f"\nReaction: {rxn.name}")
    print(f"Precedents: {len(rxn.precedents)}")

    for prec in rxn.precedents[:2]:
        print(f"  - {prec.name}")
        print(f"    Similarity: {prec.similarity:.3f}")
        print(f"    Data: {prec.data}")

print("\n=== EnzymeMap Precedents ===")
for rxn in em_reactions[:3]:
    print(f"\nReaction: {rxn.name}")
    print(f"Precedents: {len(rxn.precedents)}")

    for prec in rxn.precedents[:2]:
        print(f"  - Similarity: {prec.similarity:.3f}")
        print(f"    Reaction: {prec.rxn_smi}")
```

## Data Model Reference

### Reaction Object

```python
# Accessing reaction attributes
reaction = reactions[0]

print(reaction.product)          # Product SMILES
print(reaction.substrates)       # List of substrate SMILES
print(reaction.score)            # Expander-specific score
print(reaction.name)             # Reaction name
print(reaction.rxn_type)         # Expander type
print(reaction.rxn_domain)       # 'biocatalysis', 'chemistry', 'biosynthesis'
print(reaction.unique_id)        # UUID
print(reaction.template_metadata)  # Template info
print(reaction.precedents)       # List of Precedent objects

# Methods
rxn_smiles = reaction.reaction_smiles()  # "SUBSTRATE>>PRODUCT"
complexity = reaction.get_complexity_change()
similarity = reaction.get_similarity_score()
rxn_dict = reaction.to_dict()
```

### Pathway Object

```python
# Accessing pathway attributes
pathway = pathways[0]

print(pathway.reactions)         # List of Reaction objects
print(pathway.target_smi)        # Target molecule
print(pathway.pathway_length)    # Number of steps
print(pathway.product_smis)      # All products
print(pathway.substrate_smis)    # All substrates
print(pathway.all_smis)          # All molecules
print(pathway.tree)              # Tree structure dict

# Methods
end_molecules = pathway.end_smis()
pa_routes = pathway.get_pa_route(sme)
saved_data = pathway.save()
reaction = pathway.get_reaction_with_product('CCCC=O')
```

## Configuration and Optimization

### MCTS Configuration

```python
from rbc2 import MCTS

mcts = MCTS(target_smi, expanders)

# Time limits
mcts.config.max_search_time = 60        # Max seconds
mcts.config.max_iterations = 5000       # Max iterations

# Search depth
mcts.config.max_depth = 12              # Max pathway steps

# UCB parameters
mcts.config.c_exploration = 1.414       # Exploration vs exploitation

# Pathway filtering
mcts.config.min_pathway_score = 0.1     # Score threshold
```

### Expander Configuration

```python
from rbc2.configs import Expansion_Config

config = Expansion_Config()

# General settings
config.max_reactions = 50               # Max reactions per expansion
config.min_reaction_score = 0.2         # Score threshold
config.max_time = 120                   # Max expansion time

# Precedent settings
config.return_full_precedents = True
config.max_precedents = 10

# Apply to expander
expander = RetroBioCatExpander(config=config)
```

## Performance Optimization

### For Fast Screening

```python
# Quick screening setup
expanders = get_expanders(['retrobiocat'])  # Single expander

mcts = MCTS(target_smi, expanders)
mcts.config.max_search_time = 10     # Short time
mcts.config.max_depth = 6            # Shallow search
mcts.config.c_exploration = 2.0      # More exploration
mcts.run()
```

### For Thorough Search

```python
# Comprehensive search
expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms',
                          'retrorules', 'aizynthfinder'])

mcts = MCTS(target_smi, expanders)
mcts.config.max_search_time = 300    # 5 minutes
mcts.config.max_iterations = 10000
mcts.config.max_depth = 15
mcts.config.c_exploration = 1.0      # More exploitation
mcts.run()
```

### Parallel Processing

```python
from concurrent.futures import ProcessPoolExecutor
from rbc2 import MCTS, get_expanders

def plan_route(target_smi):
    expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
    mcts = MCTS(target_smi, expanders)
    mcts.config.max_search_time = 30
    mcts.run()
    return mcts.get_solved_pathways()

targets = ['CC(O)C(=O)O', 'CCO', 'c1ccccc1O']

with ProcessPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(plan_route, targets))

for target, pathways in zip(targets, results):
    print(f"{target}: {len(pathways)} pathways")
```

## Integration with Other Tools

### With RDKit

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rbc2 import MCTS, get_expanders

# Pre-process target
target_mol = Chem.MolFromSmiles('c1ccccc1CCO')
target_canonical = Chem.MolToSmiles(target_mol)

# Check properties
mw = Descriptors.MolWt(target_mol)
logp = Descriptors.MolLogP(target_mol)

print(f"Target MW: {mw:.2f}, LogP: {logp:.2f}")

# Run planning
expanders = get_expanders(['retrobiocat'])
mcts = MCTS(target_canonical, expanders)
mcts.run()

# Analyze products
pathways = mcts.get_solved_pathways()
for pathway in pathways:
    for smi in pathway.end_smis():
        mol = Chem.MolFromSmiles(smi)
        complexity = Descriptors.BertzCT(mol)
        print(f"  {smi}: Complexity {complexity:.1f}")
```

### With pandas for Analysis

```python
import pandas as pd
from rbc2 import MCTS, get_expanders

# Run planning
expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
mcts = MCTS(target_smi, expanders)
mcts.run()

pathways = mcts.get_all_pathways()

# Create DataFrame
data = []
for i, pathway in enumerate(pathways):
    bio_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'biocatalysis')
    chem_rxns = sum(1 for r in pathway.reactions if r.rxn_domain == 'chemistry')

    data.append({
        'pathway_id': i,
        'length': pathway.pathway_length,
        'total_reactions': len(pathway.reactions),
        'bio_reactions': bio_rxns,
        'chem_reactions': chem_rxns,
        'num_molecules': len(pathway.all_smis)
    })

df = pd.DataFrame(data)
print(df.describe())
```

## Troubleshooting

### Common Issues

**1. Installation fails on ARM Mac**
```bash
# Solution: Install HDF5 dependencies first
brew install hdf5 c-blosc
export HDF5_DIR=/opt/homebrew/opt/hdf5
export BLOSC_DIR=/opt/homebrew/opt/c-blosc
pip install tables  # Test installation
```

**2. First run downloads fail**
```python
# Expanders download data on first use
# Ensure internet connection and sufficient disk space (~500MB)

from rbc2 import RetroBioCatExpander
expander = RetroBioCatExpander()  # Downloads data automatically
```

**3. No pathways found**
```python
# Try:
# 1. Increase search time
mcts.config.max_search_time = 120

# 2. Increase max depth
mcts.config.max_depth = 15

# 3. Use more expanders
expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms', 'aizynthfinder'])

# 4. Lower similarity thresholds
expander = RetroBioCatExpander(similarity_cutoff=0.4)
```

**4. Slow performance**
```python
# Solutions:
# 1. Reduce number of expanders
expanders = get_expanders(['retrobiocat'])  # Instead of all

# 2. Reduce cutoff numbers
expander = EnzymeMapExpander(cutoff_number=20)  # Instead of 50

# 3. Set time limits
mcts.config.max_search_time = 30
```

**5. Memory issues**
```python
# Reduce pathway storage
mcts.config.max_iterations = 1000  # Lower iteration limit

# Process in batches
for target in targets[:10]:  # Process 10 at a time
    mcts = MCTS(target, expanders)
    # ... process and save
    del mcts  # Free memory
```

## Best Practices

1. **Expander Selection:**
   - Use RetroBioCat for well-known enzymatic transformations
   - Add EnzymeMap for broader enzyme coverage
   - Include AIZynthfinder for hybrid routes
   - Use RetroRules for metabolic engineering

2. **Search Configuration:**
   - Start with short search times (15-30s) for quick screening
   - Increase to 60-120s for thorough searches
   - Adjust `c_exploration` based on needs (higher = more exploration)

3. **Starting Material Evaluation:**
   - Always use appropriate SME (Commercial, E. coli, custom)
   - Create custom SMEs for lab-specific inventories
   - Consider chirality requirements

4. **Pathway Analysis:**
   - Filter pathways by domain (bio-only, hybrid)
   - Check end molecule availability
   - Analyze literature precedents for validation

5. **Integration:**
   - Use RDKit for molecule validation and property calculation
   - Export to pandas for statistical analysis
   - Save pathways as JSON for reproducibility

## Quick Reference

### Essential Imports

```python
# Core modules
from rbc2 import (
    MCTS,
    get_expanders,
    RetroBioCatExpander,
    EnzymeMapExpander,
    BKMSExpander,
    RetroRulesExpander,
    AIZynthfinderExpander,
    CommercialSME,
    EcoliSME,
    load_pathway
)

# Configuration
from rbc2.configs import Expansion_Config
```

### Quick Start Templates

**Single-step:**
```python
from rbc2 import RetroBioCatExpander

expander = RetroBioCatExpander()
reactions = expander.get_reactions('CCCC=O')
```

**Multi-step:**
```python
from rbc2 import MCTS, get_expanders

expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
mcts = MCTS(target_smi, expanders)
mcts.config.max_search_time = 30
mcts.run()
pathways = mcts.get_solved_pathways()
```

**With availability check:**
```python
from rbc2 import MCTS, get_expanders, CommercialSME

expanders = get_expanders(['retrobiocat'])
sme = CommercialSME()
mcts = MCTS(target_smi, expanders, starting_material_evaluator=sme)
mcts.run()
```

## References

**Original Paper:**
- Finnigan, W., Hepworth, L. J., Flitsch, S. L. & Turner, N. J. RetroBioCat as a computer-aided synthesis planning tool for biocatalytic reactions and cascades. Nature Catalysis 4, 98–104 (2021)

**GitHub Repository:**
- https://github.com/willfinnigan/RetroBioCat-2

**Documentation:**
- https://retrobiocat-2.readthedocs.io/

**Related Tools:**
- EnzymeMap: https://github.com/hesther/enzymemap
- RetroRules: https://retrorules.org/
- AIZynthfinder: https://github.com/MolecularAI/aizynthfinder
