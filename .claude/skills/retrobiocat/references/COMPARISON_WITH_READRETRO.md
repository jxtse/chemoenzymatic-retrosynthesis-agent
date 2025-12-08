# RetroBioCat 2.0 vs READRetro: Comparison Guide

## Quick Decision Guide

### Use READRetro when:
- ✅ Focus on **natural product biosynthesis**
- ✅ Need **KEGG pathway** integration
- ✅ Working with **secondary metabolites**
- ✅ Want **dual-view ensemble** (template-based + template-free)
- ✅ Prefer **pre-configured** setup (models already trained)
- ✅ Need **fast single-run** predictions

### Use RetroBioCat 2.0 when:
- ✅ Need **hybrid chemistry + biocatalysis** routes
- ✅ Want to **compare multiple approaches** (7 expanders)
- ✅ Require **custom expanders** or evaluators
- ✅ Need **literature precedent** searching (BRENDA, BKMS)
- ✅ Want **modular, extensible** framework
- ✅ Building **custom workflows** or integrations

### Use Both when:
- ✅ Comprehensive biosynthesis planning project
- ✅ Comparing different retrosynthesis strategies
- ✅ Research requiring multiple validation approaches
- ✅ Building integrated synthesis planning pipelines

## Feature Comparison

| Feature | READRetro | RetroBioCat 2.0 |
|---------|-----------|----------------|
| **Primary Focus** | Natural product biosynthesis | Biocatalysis + chemistry hybrid |
| **Approach** | Dual-view (Retroformer + G2S) | Multiple expanders (7 options) |
| **Data Source** | KEGG pathways | BRENDA, BKMS, RetroRules, custom |
| **Algorithm** | Retro* search | Monte Carlo Tree Search (MCTS) |
| **Setup** | Pre-trained models (~600MB) | Install + download on first use (~500MB) |
| **Installation** | conda environment | pip from GitHub |
| **Interface** | Python scripts + web | Python API only |
| **Customization** | Fixed dual-view | Highly modular |
| **Template Types** | KEGG reactions | Multiple databases |
| **Literature Precedents** | KEGG reactions | BRENDA, BKMS, custom |
| **Starting Materials** | KEGG building blocks | Commercial, E. coli, custom |
| **Output Format** | Text routes | Pathway objects with metadata |
| **Speed** | Fast (~23s per molecule) | Configurable (10s - 5min) |
| **Chemistry Support** | Enzymatic only | Enzymatic + organic chemistry |

## Use Case Examples

### Example 1: Natural Product Biosynthesis

**Scenario:** Plan biosynthesis route for a terpene natural product

```python
# READRetro approach
from readretro import run_single_molecule
pathways = run_single_molecule('CC(C)=CCCC(C)=CCO')  # Geraniol
# Uses KEGG pathways, MEP/MVA pathway knowledge

# RetroBioCat approach
from rbc2 import MCTS, get_expanders
expanders = get_expanders(['retrobiocat', 'retrorules', 'bkms'])
mcts = MCTS('CC(C)=CCCC(C)=CCO', expanders)
mcts.run()
pathways = mcts.get_solved_pathways()
# Explores multiple enzymatic databases
```

**Recommendation:** Use **READRetro** for KEGG-pathway-based routes, **RetroBioCat** for broader exploration including non-KEGG enzymes.

### Example 2: Pharmaceutical Intermediate

**Scenario:** Hybrid synthesis of an ibuprofen intermediate

```python
# READRetro: Limited to enzymatic steps
# Not ideal for hybrid routes

# RetroBioCat: Supports hybrid
from rbc2 import MCTS, get_expanders
expanders = get_expanders(['retrobiocat', 'aizynthfinder', 'askcos'])
mcts = MCTS('CC(C)Cc1ccc(cc1)C(C)C(=O)O', expanders)
mcts.run()

# Filter for hybrid routes
for pathway in mcts.get_solved_pathways():
    bio_steps = sum(1 for r in pathway.reactions if r.rxn_domain == 'biocatalysis')
    chem_steps = sum(1 for r in pathway.reactions if r.rxn_domain == 'chemistry')
    print(f"Bio: {bio_steps}, Chem: {chem_steps}")
```

**Recommendation:** Use **RetroBioCat 2.0** for hybrid chemistry + biocatalysis.

### Example 3: Metabolic Engineering

**Scenario:** Design E. coli pathway for non-natural product

```python
# READRetro: KEGG-based
from readretro import run_single_molecule
pathways = run_single_molecule('CC(O)C(=O)O')  # Lactic acid
# Will use KEGG building blocks

# RetroBioCat: E. coli specific
from rbc2 import MCTS, get_expanders, EcoliSME
expanders = get_expanders(['retrorules', 'bkms', 'enzymemap'])
sme = EcoliSME()
mcts = MCTS('CC(O)C(=O)O', expanders, starting_material_evaluator=sme)
mcts.run()
# Filters to E. coli metabolites
```

**Recommendation:** Use **RetroBioCat 2.0** with EcoliSME for organism-specific pathway design.

### Example 4: Literature Validation

**Scenario:** Find literature precedents for enzymatic reactions

```python
# READRetro: KEGG reactions as precedents
from readretro import run_single_molecule
pathways = run_single_molecule('COc1cc(ccc1O)C=O')  # Vanillin
# Shows KEGG reaction IDs

# RetroBioCat: Multiple literature sources
from rbc2 import RetroBioCatExpander, EnzymeMapExpander
rbc = RetroBioCatExpander(search_literature_precedent=True)
em = EnzymeMapExpander(enable_precedent_search=True)

rbc_reactions = rbc.get_reactions('COc1cc(ccc1O)C=O')
for rxn in rbc_reactions:
    for prec in rxn.precedents:
        print(f"{prec.name}: {prec.data}")  # BRENDA data, papers, etc.
```

**Recommendation:** Use **RetroBioCat 2.0** for detailed literature precedent analysis from BRENDA and other sources.

## Performance Comparison

### Speed

**READRetro:**
- Single molecule: ~23 seconds
- Fixed computation time
- Optimized for throughput

**RetroBioCat 2.0:**
- Configurable: 10s (quick) to 300s (thorough)
- Depends on expanders used
- Flexible time/quality tradeoff

### Memory Usage

**READRetro:**
- ~2-4GB RAM for models
- Consistent usage

**RetroBioCat 2.0:**
- ~1-3GB depending on expanders
- More with multiple expanders simultaneously

### Scalability

**READRetro:**
- Batch mode for multiple targets
- Multiprocessing support built-in

**RetroBioCat 2.0:**
- MCTS naturally parallelizable
- Manual parallel processing implementation needed

## Integration Possibilities

### Using Both Tools Together

```python
# Strategy: Use READRetro for KEGG-based initial screen,
#          then RetroBioCat for detailed exploration

# Step 1: Quick KEGG-based screening with READRetro
from readretro import run_single_molecule
readretro_paths = run_single_molecule(target_smi)

# Step 2: Detailed exploration with RetroBioCat
from rbc2 import MCTS, get_expanders
expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms'])
mcts = MCTS(target_smi, expanders)
mcts.config.max_search_time = 60
mcts.run()
rbc_paths = mcts.get_solved_pathways()

# Step 3: Compare and validate
# - KEGG pathway coverage (READRetro)
# - Literature precedents (RetroBioCat)
# - Route diversity (both)
```

### Complementary Workflow

1. **Initial Screening (READRetro)**
   - Fast KEGG-based route generation
   - Identify feasible pathway classes
   - Benchmark against known biosynthetic routes

2. **Detailed Planning (RetroBioCat 2.0)**
   - Explore non-KEGG enzymatic options
   - Add hybrid chemistry steps if needed
   - Search detailed literature precedents
   - Validate with commercial availability

3. **Validation**
   - Compare routes from both tools
   - Cross-reference KEGG pathways with BRENDA data
   - Select optimal route based on:
     - Pathway length
     - Literature precedents
     - Commercial availability
     - Organism compatibility

## Installation and Setup

### READRetro Setup
```bash
# 1. Create environment
conda create -n readretro python=3.8 -y
conda activate readretro

# 2. Install dependencies
pip install torch==1.12.0+cu113 --extra-index-url https://download.pytorch.org/whl/cu113
pip install rdkit-pypi easydict pandas tqdm numpy==1.22 OpenNMT-py==2.3.0 networkx==2.5 scipy gin-config==0.4.0

# 3. Download and setup data (from Zenodo)
# Copy models and data to correct locations
```

### RetroBioCat 2.0 Setup
```bash
# 1. Create environment
conda create -n retrobiocat python=3.9 -y
conda activate retrobiocat

# 2. Install (downloads dependencies automatically)
pip install git+https://github.com/willfinnigan/RetroBioCat-2.git

# 3. First run downloads required data automatically
python -c "from rbc2 import RetroBioCatExpander; RetroBioCatExpander()"
```

## Which Should You Learn First?

### For Beginners:
**Start with READRetro** if:
- You're working specifically with natural products
- You want a straightforward, pre-configured setup
- You need results quickly with minimal configuration

**Start with RetroBioCat 2.0** if:
- You need flexibility and customization
- You're comfortable with Python APIs
- You want to understand modular retrosynthesis frameworks

### For Advanced Users:
**Use both** to leverage:
- READRetro's optimized KEGG pathway knowledge
- RetroBioCat's extensible architecture
- Complementary data sources (KEGG vs BRENDA/BKMS)
- Different search algorithms (Retro* vs MCTS)

## Summary Table

| Criterion | READRetro | RetroBioCat 2.0 | Winner |
|-----------|-----------|----------------|--------|
| **Natural Products** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | READRetro |
| **Hybrid Routes** | ⭐ | ⭐⭐⭐⭐⭐ | RetroBioCat |
| **Ease of Use** | ⭐⭐⭐⭐ | ⭐⭐⭐ | READRetro |
| **Customization** | ⭐⭐ | ⭐⭐⭐⭐⭐ | RetroBioCat |
| **Speed** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | READRetro |
| **Data Coverage** | ⭐⭐⭐⭐ (KEGG) | ⭐⭐⭐⭐⭐ (Multiple) | RetroBioCat |
| **Precedents** | ⭐⭐⭐ (KEGG) | ⭐⭐⭐⭐⭐ (BRENDA+) | RetroBioCat |
| **Setup** | ⭐⭐⭐ | ⭐⭐⭐⭐ | RetroBioCat |
| **Documentation** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Tie |
| **Community** | ⭐⭐⭐ | ⭐⭐⭐⭐ | RetroBioCat |

## Recommendations by Field

### Natural Product Chemistry
**Primary:** READRetro
**Secondary:** RetroBioCat (for non-KEGG enzymes)

### Pharmaceutical Development
**Primary:** RetroBioCat 2.0 (hybrid routes)
**Secondary:** READRetro (bio-only validation)

### Metabolic Engineering
**Primary:** RetroBioCat 2.0 (with EcoliSME)
**Secondary:** READRetro (KEGG pathway reference)

### Synthetic Biology
**Primary:** Both (complementary)
**Secondary:** Custom integration

### Academic Research
**Primary:** Both (comprehensive validation)
**Secondary:** Method comparison studies

## Final Thoughts

Both tools are excellent for biosynthesis planning but serve different needs:

- **READRetro** excels at rapid, KEGG-informed natural product biosynthesis with proven pathway knowledge
- **RetroBioCat 2.0** shines in flexible, modular planning with hybrid routes and extensive literature coverage

For comprehensive biosynthesis projects, consider using both tools to leverage their complementary strengths.
