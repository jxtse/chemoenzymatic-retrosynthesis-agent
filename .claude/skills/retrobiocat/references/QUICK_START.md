# RetroBioCat 2.0 Quick Start Guide

## 5-Minute Quick Start

### Installation
```bash
conda create -n retrobiocat python=3.9 -y
conda activate retrobiocat
pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### Your First Retrosynthesis

```python
from rbc2 import MCTS, get_expanders

# Define target molecule
target = 'CC(O)C(=O)O'  # Lactic acid

# Get expanders
expanders = get_expanders(['retrobiocat', 'aizynthfinder'])

# Run MCTS
mcts = MCTS(target, expanders)
mcts.config.max_search_time = 30
mcts.run()

# Get results
pathways = mcts.get_solved_pathways()
print(f"Found {len(pathways)} pathways")

# Show best pathway
if pathways:
    best = pathways[0]
    print(f"\nBest pathway: {best.pathway_length} steps")
    for i, rxn in enumerate(best.reactions, 1):
        print(f"{i}. {rxn.reaction_smiles()}")
```

## Common Tasks

### Task 1: Single-Step Retrosynthesis
```python
from rbc2 import RetroBioCatExpander

expander = RetroBioCatExpander()
reactions = expander.get_reactions('CCCC=O')

for rxn in reactions[:5]:
    print(f"{rxn.reaction_smiles()} (score: {rxn.score:.2f})")
```

### Task 2: Compare Multiple Expanders
```python
from rbc2 import get_expanders

target = 'c1ccc(cc1)CCO'
expanders = get_expanders(['retrobiocat', 'enzymemap', 'aizynthfinder'])

for exp in expanders:
    rxns = exp.get_reactions(target)
    print(f"{exp.__class__.__name__}: {len(rxns)} reactions")
```

### Task 3: Check Commercial Availability
```python
from rbc2 import CommercialSME

sme = CommercialSME()
available, info = sme.eval('CCO')
print(f"Available: {available}")
```

### Task 4: Bio-Only Routes
```python
from rbc2 import MCTS, get_expanders

expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms'])
mcts = MCTS(target, expanders)
mcts.run()

# Filter bio-only
bio_only = [p for p in mcts.get_solved_pathways()
            if all(r.rxn_domain == 'biocatalysis' for r in p.reactions)]
```

## Configuration Presets

### Quick Screening
```python
mcts.config.max_search_time = 10
mcts.config.max_depth = 6
expanders = get_expanders(['retrobiocat'])
```

### Balanced
```python
mcts.config.max_search_time = 30
mcts.config.max_depth = 10
expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
```

### Thorough
```python
mcts.config.max_search_time = 120
mcts.config.max_depth = 15
expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms', 'aizynthfinder'])
```

## Helper Script Usage

```bash
# Single-step
python scripts/retrobiocat_helper.py single -s "CCCC=O" -e retrobiocat

# Multi-step
python scripts/retrobiocat_helper.py multi -s "CC(O)C(=O)O" -t 60

# Check availability
python scripts/retrobiocat_helper.py check -s "CCO" "c1ccccc1O"

# Batch processing
python scripts/retrobiocat_helper.py batch -f targets.txt -o results.json
```

## Troubleshooting

**Problem:** Installation fails on ARM Mac
```bash
brew install hdf5 c-blosc
export HDF5_DIR=/opt/homebrew/opt/hdf5
export BLOSC_DIR=/opt/homebrew/opt/c-blosc
pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

**Problem:** No pathways found
```python
# Increase time and depth
mcts.config.max_search_time = 120
mcts.config.max_depth = 15

# Use more expanders
expanders = get_expanders(['retrobiocat', 'enzymemap', 'aizynthfinder'])
```

**Problem:** Slow performance
```python
# Use fewer expanders
expanders = get_expanders(['retrobiocat'])

# Reduce cutoff
expander = EnzymeMapExpander(cutoff_number=20)
```

## Next Steps

1. Read full SKILL.md for comprehensive examples
2. Explore different expander combinations
3. Customize starting material evaluators
4. Integrate with your workflow

## Quick Reference

**Import essentials:**
```python
from rbc2 import (
    MCTS, get_expanders,
    RetroBioCatExpander, EnzymeMapExpander,
    CommercialSME, EcoliSME
)
```

**Basic workflow:**
1. Create expanders
2. Initialize MCTS with target
3. Configure search
4. Run and get pathways
5. Analyze results

**Documentation:**
- Full skill: `.claude/skills/retrobiocat/SKILL.md`
- References: `.claude/skills/retrobiocat/references/`
- Online: https://retrobiocat-2.readthedocs.io/
