---
name: readretro
description: "Natural product biosynthesis planning with retrieval-augmented dual-view retrosynthesis. Plan enzymatic synthesis routes from target molecules, multi-step retrosynthesis, KEGG pathway integration, for natural product and biocatalysis research."
---

# READRetro: Natural Product Biosynthesis Planning

## Overview

READRetro (Retrieval-Augmented Dual-View Retrosynthesis) is a specialized retrosynthesis tool designed for natural product biosynthesis planning. It combines template-based (Retroformer) and template-free (Graph2SMILES) models with KEGG pathway knowledge to predict biologically feasible enzymatic synthesis routes.

**Core Capabilities:**
- Multi-step retrosynthesis planning for natural products
- Dual-view prediction (template-based + template-free)
- KEGG pathway knowledge integration
- Biochemical building block recognition
- Enzymatic reaction feasibility assessment
- Compatible with both GPU and CPU modes

**Key Distinction:** READRetro is optimized for **biosynthesis routes** (enzymatic reactions), not general chemical synthesis. Use for natural products, secondary metabolites, and biocatalysis applications.

## When to Use This Skill

This skill should be used when:

- "Plan biosynthesis pathway for this molecule"
- "Predict enzymatic synthesis route"
- "Natural product retrosynthesis"
- "Find biosynthetic pathway using KEGG"
- "Multi-step enzymatic synthesis planning"
- Tasks involving natural products or secondary metabolites
- Biocatalysis and metabolic engineering applications
- SMILES-based biosynthesis route prediction

## Installation and Environment Setup

### Environment Status

The `readretro` conda environment has been pre-configured with all dependencies.

**Environment Details:**
- Python 3.8
- PyTorch 1.12.0+cu113 (CUDA 11.3)
- RDKit (via rdkit-pypi)
- OpenNMT-py 2.3.0
- Other dependencies: pandas, numpy, networkx, scipy, gin-config

### Activate Environment

```bash
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
```

### Verify Installation

```bash
python -c "import torch; import rdkit; print('PyTorch:', torch.__version__); print('CUDA:', torch.cuda.is_available())"
```

## Data and Model Files

**Required Files (Already Configured):**
- `READRetro/retroformer/saved_models/biochem.pt` - Retroformer model
- `READRetro/g2s/saved_models/biochem.pt` - Graph2SMILES model
- `READRetro/data/` - KEGG reaction database and building blocks
  - `kegg_reaction.pickle` - KEGG reaction templates
  - `kegg_neutral_iso_smi.csv` - KEGG compound structures
  - `building_block.csv` - Biosynthetic building blocks
  - `pathways.pickle` - KEGG pathway information

## Core Workflows

### Workflow 1: Single Molecule Retrosynthesis

**Use Case:** Plan biosynthesis route for one target molecule

**Input:** SMILES string of target molecule

**Command:**
```bash
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro

# Run with GPU
CUDA_VISIBLE_DEVICES=0 python run.py 'O=C1C=C2C=CC(O)CC2O1'

# Run with CPU (slower)
python run.py 'O=C1C=C2C=CC(O)CC2O1'
```

**Output:**
- Retrosynthesis tree printed to console
- Route saved to `result/` directory (if configured)
- Shows predicted precursors and enzymatic reactions

**Example:**
```bash
# Natural product example (paracetamol derivative)
python run.py 'CC(=O)Nc1ccc(O)cc1'

# Terpene example
python run.py 'CC(C)=CCCC(C)=C'

# Flavonoid example
python run.py 'O=C1C=C(c2ccccc2)Oc2cc(O)cc(O)c12'
```

### Workflow 2: Batch Multi-Step Planning

**Use Case:** Plan routes for multiple molecules in parallel

**Input:** Molecule list defined in `run_mp.py` or data files

**Command:**
```bash
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro

# Multi-threaded batch processing
CUDA_VISIBLE_DEVICES=0 python run_mp.py
```

**Configuration (in run_mp.py):**
```python
num_threads = 4  # Adjust based on GPU memory
beam_size = 10   # Number of routes to explore
max_depth = 15   # Maximum retrosynthesis steps
```

**Output:**
- Results saved to `result/` directory
- One file per molecule with predicted routes
- Includes reaction templates and confidence scores

### Workflow 3: Single-Step Retrosynthesis Evaluation

**Use Case:** Evaluate one-step retrosynthesis accuracy

**Commands:**
```bash
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro

# Ensemble model (both Retroformer + Graph2SMILES)
CUDA_VISIBLE_DEVICES=0 python eval_single.py

# Retroformer only
CUDA_VISIBLE_DEVICES=0 python eval_single.py -m retroformer

# Graph2SMILES only (with 200 samples)
CUDA_VISIBLE_DEVICES=0 python eval_single.py -m g2s -s 200
```

**Output:**
- Top-k accuracy metrics (k=1,3,5,10)
- Per-molecule prediction results
- Performance on test set

### Workflow 4: Route Evaluation

**Use Case:** Evaluate predicted multi-step routes

**Command:**
```bash
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro

# Evaluate routes from batch planning
python eval.py result/output_file.txt
```

**Metrics:**
- Success rate (routes reaching building blocks)
- Average route length
- KEGG pathway coverage

## Key Parameters and Configuration

### Search Parameters

**In `run.py` and `run_mp.py`:**

```python
# Retrosynthesis search parameters
beam_size = 10          # Number of parallel routes (default: 10)
max_depth = 15          # Maximum synthesis steps (default: 15)
num_path = 100          # Total paths to explore (default: 100)

# Model parameters
retroformer_topk = 50   # Top-k predictions from Retroformer
g2s_topk = 50          # Top-k predictions from Graph2SMILES

# Building blocks
use_kegg_building_blocks = True  # Use KEGG starting materials
```

### Model Selection

**Available Models:**
- `retroformer` - Template-based (faster, more interpretable)
- `g2s` - Template-free (more flexible, handles novel reactions)
- `ensemble` - Both models (default, best performance)

**Model Checkpoints:**
- Biochemical dataset: `biochem.pt`
- Clean dataset: `clean.pt` (USPTO-based)

## Common Use Cases and Examples

### Use Case 1: Natural Product Biosynthesis

```bash
# Target: Resveratrol (stilbenoid natural product)
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
conda activate readretro
python run.py 'C=1C=C(C=CC1O)C=CC2=CC(=CC(=C2)O)O'
```

**Expected Output:**
- Phenylpropanoid pathway routes
- Phenylalanine as starting material
- Enzymatic steps (PAL, 4CL, STS)

### Use Case 2: Terpene Synthesis

```bash
# Target: Geraniol (monoterpene)
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
conda activate readretro
python run.py 'CC(C)=CCCC(C)=CCO'
```

**Expected Output:**
- Mevalonate or MEP pathway
- Isoprene units (IPP, DMAPP)
- Terpene synthase reactions

### Use Case 3: Screening Compound Library

```python
# Create input file: compounds.txt
# One SMILES per line

# Modify run_mp.py to read from compounds.txt
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
conda activate readretro
python run_mp.py

# Results saved to result/ directory
```

### Use Case 4: Compare Routes from Different Models

```bash
# Retroformer route
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro
conda activate readretro
python run.py 'YOUR_SMILES' --model retroformer > routes_rf.txt

# Graph2SMILES route
python run.py 'YOUR_SMILES' --model g2s > routes_g2s.txt

# Compare route diversity and feasibility
```

## Performance Optimization

### GPU Memory Management

```bash
# Reduce beam size for limited GPU memory
python run.py 'SMILES' --beam_size 5

# Use CPU if GPU unavailable
python run.py 'SMILES'  # Automatically uses CPU
```

### Batch Processing

```python
# In run_mp.py, adjust threads based on GPU:
# GPU with 12GB VRAM: num_threads = 4
# GPU with 24GB VRAM: num_threads = 8
# CPU: num_threads = 2
```

### Speed vs Accuracy Trade-offs

```bash
# Fast (lower accuracy): Small beam, few paths
python run.py 'SMILES' --beam_size 5 --num_path 50

# Balanced (default)
python run.py 'SMILES' --beam_size 10 --num_path 100

# Thorough (slower, higher accuracy)
python run.py 'SMILES' --beam_size 20 --num_path 200
```

## Interpreting Results

### Route Output Format

```
Product: [SMILES]
├── Reaction 1: [Template ID]
│   ├── Reactant 1: [SMILES]
│   └── Reactant 2: [SMILES]
├── Reaction 2: [Template ID]
│   └── Reactant: [SMILES]
...
└── Building Blocks: [List of KEGG compounds]
```

### Confidence Scores

- **Template Score**: Retroformer confidence (0-1)
- **Path Score**: Cumulative route confidence
- **KEGG Coverage**: Fraction of reactions matching KEGG

### Success Criteria

A route is considered successful if:
1. All leaf nodes are KEGG building blocks
2. All reactions have valid templates
3. Route length ≤ max_depth

## Troubleshooting

### Common Issues

**1. Import Error (scipy, gin, etc.)**
```bash
conda activate readretro
pip install scipy gin-config==0.4.0
```

**2. CUDA Out of Memory**
```bash
# Reduce beam size or number of threads
python run.py 'SMILES' --beam_size 5
# Or use CPU
python run.py 'SMILES'
```

**3. Model Files Not Found**
```bash
# Verify model paths
ls READRetro/retroformer/saved_models/biochem.pt
ls READRetro/g2s/saved_models/biochem.pt
```

**4. Invalid SMILES**
```python
# Pre-validate SMILES with RDKit
from rdkit import Chem
mol = Chem.MolFromSmiles('YOUR_SMILES')
if mol is None:
    print("Invalid SMILES")
else:
    canonical_smiles = Chem.MolToSmiles(mol)
    # Use canonical_smiles
```

**5. Slow Performance**
- Ensure CUDA is available: `python -c "import torch; print(torch.cuda.is_available())"`
- Reduce `num_path` and `beam_size`
- Use single model instead of ensemble

## Integration with Other Tools

### With RDKit (Molecule Validation)

```python
from rdkit import Chem
import subprocess

def run_readretro(smiles):
    # Validate and canonicalize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"

    canonical = Chem.MolToSmiles(mol)

    # Run READRetro
    cmd = f"cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro && conda activate readretro && python run.py '{canonical}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout
```

### With KEGG Database

READRetro internally uses KEGG data (`kegg_reaction.pickle`, `kegg_neutral_iso_smi.csv`). To update:

```python
# Access KEGG reaction templates
import pickle
with open('data/kegg_reaction.pickle', 'rb') as f:
    reactions = pickle.load(f)

# Access KEGG compounds
import pandas as pd
kegg_compounds = pd.read_csv('data/kegg_neutral_iso_smi.csv')
```

### With Pathway Databases

```python
# Load pathway information
import pickle
with open('data/pathways.pickle', 'rb') as f:
    pathways = pickle.load(f)

# Map predicted reactions to pathways
for reaction in predicted_route:
    pathway_ids = pathways.get(reaction['template_id'], [])
    print(f"Pathway: {pathway_ids}")
```

## Best Practices

1. **Input Preparation:**
   - Canonicalize SMILES with RDKit before input
   - Remove salts and neutralize charges
   - Verify molecular structure validity

2. **Parameter Selection:**
   - Start with default parameters (beam=10, paths=100)
   - Increase for complex molecules or higher confidence
   - Decrease for quick screening or limited resources

3. **Result Validation:**
   - Check KEGG coverage of predicted routes
   - Verify biochemical feasibility of enzymes
   - Cross-reference with literature pathways

4. **Batch Processing:**
   - Group similar molecules for efficient processing
   - Monitor GPU memory usage
   - Save intermediate results

5. **Reproducibility:**
   - Record random seeds if using stochastic components
   - Document model versions and parameters
   - Save complete route outputs

## References and Resources

**Paper:**
- Lee, S., et al. (2023). "READRetro: Natural Product Biosynthesis Planning with Retrieval-Augmented Dual-View Retrosynthesis." bioRxiv. https://www.biorxiv.org/content/10.1101/2023.03.21.533616v1

**Web Version:**
- https://readretro.net

**GitHub Repository:**
- https://github.com/SeulLee05/READRetro

**Data Source:**
- Zenodo: https://zenodo.org/records/11485641

**Related Tools:**
- Retroformer: Template-based retrosynthesis
- Graph2SMILES: Template-free retrosynthesis
- KEGG PATHWAY: Metabolic pathway database

## Quick Reference

### Essential Commands

```bash
# Activate environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate readretro
cd /home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro

# Single molecule
python run.py 'SMILES_STRING'

# Batch processing
python run_mp.py

# Evaluation
python eval_single.py
python eval.py result/output.txt
```

### Directory Structure

```
READRetro/
├── run.py                  # Single molecule planning
├── run_mp.py              # Batch processing
├── eval_single.py         # Single-step evaluation
├── eval.py                # Multi-step evaluation
├── data/                  # KEGG data and building blocks
├── retroformer/           # Retroformer model
│   └── saved_models/
├── g2s/                   # Graph2SMILES model
│   └── saved_models/
└── result/                # Output directory
```

### Model Files

- Retroformer (biochem): `retroformer/saved_models/biochem.pt` (~500 MB)
- Graph2SMILES (biochem): `g2s/saved_models/biochem.pt` (~87 MB)
- KEGG reactions: `data/kegg_reaction.pickle`
- Building blocks: `data/building_block.csv`

## Environment Details

**Location:** `/home/xiejinxiang/.conda/envs/readretro`

**Key Dependencies:**
- Python: 3.8.20
- PyTorch: 1.12.0+cu113
- CUDA: 11.3
- RDKit: 2022.9.5 (via rdkit-pypi)
- OpenNMT-py: 2.3.0
- pandas: 2.0.3
- numpy: 1.22.0
- networkx: 2.5
- scipy: 1.10.1
- gin-config: 0.4.0

**Working Directory:** `/home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro`
