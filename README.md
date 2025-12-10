# Chemoenzymatic Retrosynthesis Agent

AI-powered retrosynthesis planning system combining enzymatic and chemical reactions for sustainable synthesis route design.

[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

This project provides a production-ready AI agent for chemoenzymatic retrosynthesis planning, featuring:

- **Hybrid Route Planning**: Combines biocatalytic and chemical transformations
- **Knowledge Base**: 42K+ curated enzymatic reactions and kinetic parameters
- **Intelligent Evaluation**: Automated feasibility scoring (0-10 scale)
- **Natural Language Output**: Results explained in plain language
- **Multiple Interfaces**: CLI, Python API, and batch processing

### Key Capabilities

```
Target Molecule → SMILES Conversion → RetroBioCat MCTS Planning
                                    ↓
              Commercial Availability + Enzyme Data Query
                                    ↓
              Feasibility Scoring + Plain Language Explanation
```

## Quick Start

### Prerequisites

- Python 3.10 (required, not 3.11+)
- `uv` package manager (recommended)

### Installation

```bash
# Clone repository
git clone <your-repo-url>
cd chemoenzymatic-retrosynthesis-agent

# Install dependencies
uv sync

# Install RetroBioCat 2.0 (from GitHub)
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git

# Build knowledge base (first-time only)
uv run python -m knowledge_base.cli build --config kb_config.yaml
```

### Configuration

Create a `.env` file:

```bash
cp .env.example .env
```

Edit `.env` with your API key:

```env
# Option 1: OpenRouter (Recommended - supports multiple models)
OPENROUTER_API_KEY=sk-or-v1-xxxxx
OPENROUTER_MODEL=google/gemini-flash-1.5

# Option 2: OpenAI
OPENAI_API_KEY=sk-xxxxx
OPENAI_MODEL=gpt-4o-mini

# Option 3: Azure OpenAI
AZURE_OPENAI_API_KEY=xxxxx
AZURE_OPENAI_ENDPOINT=https://xxx.openai.azure.com/
AZURE_OPENAI_DEPLOYMENT=gpt-4
```

### Usage

#### Interactive Mode

```bash
uv run python run_agent.py
```

Then enter target molecule:
```
请输入目标分子: ibuprofen
```

#### Command Line Mode

```bash
# Full planning
uv run python run_agent.py --target "ibuprofen"

# Quick check
uv run python run_agent.py --target "L-DOPA" --quick

# Pure biocatalytic route
uv run python run_agent.py --target "phenylethylamine" --no-chemistry

# Custom max steps
uv run python run_agent.py --target "aspirin" --max-steps 5
```

#### Batch Processing

```bash
# Create target list
echo -e "ibuprofen\nL-DOPA\naspirin" > targets.txt

# Batch process
uv run python run_agent.py --batch targets.txt --output results.json
```

#### Python API

```python
from agents.production_agent import RetrosynthesisAgent

# Configure LLM
llm_config = {
    "config_list": [{
        "model": "google/gemini-flash-1.5",
        "api_key": "sk-or-v1-xxxxx",
        "base_url": "https://openrouter.ai/api/v1",
    }],
    "temperature": 0.7,
}

# Create agent
agent = RetrosynthesisAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config=llm_config
)

# Plan route
result = agent.plan(target="ibuprofen", max_steps=6)
print(result)

# Quick feasibility check
quick_result = agent.quick_check("L-DOPA")
print(quick_result)
```

## Project Structure

```
chemoenzymatic-retrosynthesis-agent/
├── agents/
│   ├── production_agent.py       # Main agent implementation
│   ├── retrobiocat_tools.py      # RetroBioCat integration
│   ├── kb_tools.py               # Knowledge base query tools
│   ├── utils.py                  # Helper functions
│   └── logging_utils.py          # Session logging utilities
│
├── knowledge_base/               # Knowledge base system
│   ├── api.py                    # Query API
│   ├── builder.py                # Builder
│   ├── cli.py                    # Command-line interface
│   └── connectors/               # Database connectors (9 sources)
│
├── datasets/                     # Data sources
│   ├── Reactions_BKMS.csv        # 11MB, 42K reactions
│   ├── brenda_kcat_v3.parquet    # 3.4MB, kinetic data
│   └── EnzyExtractDB_176463.parquet  # 9.7MB, enzyme data
│
├── knowledge_base_output/        # Built knowledge base
│   └── knowledge_base.jsonl      # 76MB, 42K+ records
│
├── examples/
│   └── production_example.py     # Usage examples
│
├── run_agent.py                  # Main entry point
├── pyproject.toml                # Project dependencies
├── .env.example                  # Environment template
└── kb_config.yaml                # Knowledge base config
```

## Features

### 1. Retrosynthesis Planning

- **Single-Step Reactions**: Query enzymatic transformations for target molecule
- **Multi-Step MCTS Planning**: Hybrid biocatalytic + chemical route design
- **Commercial Availability**: Check starting material purchasability (100K+ compounds)
- **Enzyme Database**: Query EC numbers, kinetic parameters, sequences

### 2. Knowledge Base

| Database | Type | Records | Status |
|----------|------|---------|--------|
| BKMS | Local CSV | 42,539 | ✅ |
| BRENDA | Local Parquet | ~162 | ✅ |
| EnzyExtract | Local Parquet | ~170K | ✅ |
| KEGG | REST API | - | ⚠️ Optional |
| UniProt | REST API | - | ⚠️ Optional |
| PubChem | REST API | - | ⚠️ Optional |

**Total**: 42,701 unified records (76MB JSONL)

### 3. Feasibility Scoring

Routes are scored 0-10 based on:

- **Steps** (3 pts): Fewer steps score higher
- **Starting Material Availability** (3 pts): Commercial availability
- **Enzyme Availability** (2 pts): Common vs rare enzymes
- **Literature Precedent** (2 pts): Prior successful examples

**Score Interpretation**:
- **8-10**: 💚 Highly recommended
- **5-7**: 💛 Moderate feasibility, note challenges
- **1-4**: ❤️ Not recommended, consider alternatives

### 4. Plain Language Output

Example output format:

```
🛣️ Route Overview
Total 5 steps: 2 enzymatic, 3 chemical

Step 1: [Plain language explanation]
Step 2: [Plain language explanation]
...

🧬 Key Enzymes
- Enzyme Name (EC X.X.X.X): Source, efficiency, availability
...

🛒 Starting Materials
✅ Material A: Sigma-Aldrich, ~$50/100g
⚠️ Material B: Needs synthesis
...

📊 Feasibility Score: 8/10
💚 Highly recommended!

Breakdown:
- Steps (5): +2 pts
- Materials: +3 pts (all commercial)
- Enzymes: +3 pts (common)
- Precedent: +0 pts

💡 Implementation Suggestions
- Challenge: [specific issue]
- Alternative: [optional approach]
- Notes: [key points]
```

## Core Dependencies

### AI & Agent Framework
- **pyautogen** (0.2.x): Multi-agent orchestration
- **openai** (≥1.0.0): OpenAI/Azure/OpenRouter client

### Cheminformatics
- **rdkit** (≥2023.9.1): Chemical structure processing
- **pubchempy** (≥1.0.4): PubChem API
- **rbc2**: RetroBioCat 2.0 (install from GitHub)

### Data Processing
- **numpy** (1.26.x), **pandas** (≥2.2.3)
- **scipy** (≥1.11.0), **pyarrow** (≥17.0.0)

### Deep Learning
- **torch** (≥2.0.0), **tensorflow** (2.8.x)

See [DEPENDENCIES.md](DEPENDENCIES.md) for complete list and [INSTALL.md](INSTALL.md) for detailed installation.

## Examples

### Example 1: Ibuprofen Synthesis

```bash
uv run python run_agent.py --target "ibuprofen"
```

### Example 2: L-DOPA with Python API

```python
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(kb_path="knowledge_base_output/knowledge_base.jsonl")
result = agent.quick_check("L-DOPA")
print(result)
```

### Example 3: Batch Screening

```python
from agents.production_agent import RetrosynthesisAgent

targets = ["aspirin", "ibuprofen", "L-DOPA", "vanillin"]
agent = RetrosynthesisAgent(kb_path="...", llm_config={...})

for target in targets:
    result = agent.quick_check(target)
    print(f"{target}: {result}\n")
```

## Troubleshooting

### Issue: `No module named 'rbc2'`

```bash
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### Issue: Knowledge base not found

```bash
uv run python -m knowledge_base.cli build --config kb_config.yaml
```

### Issue: API key not set

```bash
# Check .env file
cat .env

# Ensure one of these is set:
# OPENROUTER_API_KEY
# OPENAI_API_KEY
# AZURE_OPENAI_API_KEY
```

### Issue: RDKit import error

```bash
uv add rdkit
```

## Advanced Usage

### Rebuild Knowledge Base

```bash
# Edit configuration
vim kb_config.yaml

# Build
uv run python -m knowledge_base.cli build --config kb_config.yaml

# Query statistics
uv run python -m knowledge_base.cli query --kb knowledge_base_output/knowledge_base.jsonl --stats
```

### Custom System Prompt

Edit `agents/production_agent.py` `_get_system_message()` method to customize agent behavior.

### Add New Tools

Add functions in `agents/utils.py`, then register in `production_agent.py`:

```python
def my_custom_tool(arg: str) -> str:
    """Custom tool description"""
    return json.dumps({"result": "..."})

# In _register_tools():
tools["my_custom_tool"] = my_custom_tool
```

## Performance

| Operation | Time |
|-----------|------|
| Quick check | 10-30s |
| Full planning (5 steps) | 1-2 min |
| Full planning (10 steps) | 2-5 min |
| Batch (10 molecules) | 10-20 min |

## Limitations

### Data Sources
- **BRENDA**: Academic use only
- **RetroBioCat**: Follow original license
- **PubChem**: Respect usage terms

### Accuracy
- Scoring is **heuristic**, not experimentally validated
- Route feasibility requires **expert chemical judgment**
- Cost estimates **not implemented**

### Performance
- Complex molecules (>10 steps) may be slow
- MCTS search time configurable (default 30s)
- Batch processing requires sufficient API quota

## Contributing

Contributions welcome! Priority areas:
- Improved feasibility scoring algorithms
- Additional database integrations
- Cost estimation functionality
- Web interface

## License

MIT License - see [LICENSE](LICENSE) for details

**Important Notes**:
- BRENDA data: Academic use only
- RetroBioCat: Follows original license
- LLM APIs: Follow provider terms

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{chemoenzymatic_retrosynthesis_agent,
  title = {Chemoenzymatic Retrosynthesis Agent},
  author = {[Your Name]},
  year = {2025},
  url = {https://github.com/[your-repo]}
}
```

## Acknowledgments

- **RetroBioCat 2.0**: W. Finnigan et al.
- **BRENDA**: Cologne University Bioinformatics Center
- **PubChem**: NCBI
- **AutoGen**: Microsoft Research

---

**Version**: 1.0.0 (Production)
**Last Updated**: 2025-12-10
**Python**: 3.10
**Package Manager**: uv
