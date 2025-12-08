# RetroBioCat 2.0 Skill References

## Primary References

### RetroBioCat Original Paper
Finnigan, W., Hepworth, L. J., Flitsch, S. L. & Turner, N. J. (2021). **RetroBioCat as a computer-aided synthesis planning tool for biocatalytic reactions and cascades**. *Nature Catalysis*, 4, 98–104.

https://doi.org/10.1038/s41929-020-00556-z

## Official Resources

- **GitHub Repository**: https://github.com/willfinnigan/RetroBioCat-2
- **Documentation**: https://retrobiocat-2.readthedocs.io/
- **Original RetroBioCat (v1)**: https://github.com/willfinnigan/RetroBioCat

## Expander-Specific References

### EnzymeMap
Heid, E., Probst, D., Green, W. H. & Madsen, G. K. H. (2023). **EnzymeMap: curation, validation and data-driven prediction of enzymatic reactions**. *Chemical Science*, 14, 14229–14242.

https://doi.org/10.1039/D3SC04546H

**GitHub**: https://github.com/hesther/enzymemap

### BKMS Expander
Levin, I., Liu, M., Voigt, C. A. & Coley, C. W. (2022). **Merging enzymatic and synthetic chemistry with computational synthesis planning**. *Nature Communications*, 13, 7747.

https://doi.org/10.1038/s41467-022-35422-y

**Data Source**: BKMS Database (Biochemical Knowledge Management System)

### RetroRules
Koch, M., Duigou, T. & Faulon, J.-L. (2020). **Reinforcement Learning for Bioretrosynthesis**. *ACS Synthetic Biology*, 9, 157–168.

https://doi.org/10.1021/acssynbio.9b00447

**Website**: https://retrorules.org/
**Database**: https://doi.org/10.5281/zenodo.5827427

### AIZynthfinder
Genheden, S., et al. (2020). **AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning**. *Journal of Cheminformatics*, 12, 70.

https://doi.org/10.1186/s13321-020-00472-1

**GitHub**: https://github.com/MolecularAI/aizynthfinder

### RingBreaker
Thakkar, A., Selmi, N., Reymond, J.-L., Engkvist, O. & Bjerrum, E. J. (2020). **"Ring Breaker": Neural Network Driven Synthesis Prediction of the Ring System Chemical Space**. *Journal of Medicinal Chemistry*, 63, 8791–8808.

https://doi.org/10.1021/acs.jmedchem.0c00761

### AskCos
Coley, C. W., et al. (2019). **A robotic platform for flow synthesis of organic compounds informed by AI planning**. *Science*, 365, eaax1566.

https://doi.org/10.1126/science.aax1566

**Website**: https://askcos.mit.edu/

## Installation Requirements

### System Requirements

**Minimum:**
- Python 3.9+
- 4GB RAM
- 2GB disk space

**Recommended:**
- Python 3.9 or 3.10
- 8GB+ RAM
- 5GB disk space (for all expander data)

### Python Dependencies

Core dependencies (installed automatically):
```
rdkit
pandas
numpy
tables (pytables)
requests
tqdm
```

Expander-specific dependencies (downloaded on first use):
- Template databases (~500MB)
- Model weights (~200MB)
- Starting material databases (~100MB)

### Platform-Specific Notes

**Linux/Ubuntu**: Works out of the box with pip installation

**macOS (Intel)**: Standard pip installation

**macOS (ARM/M1/M2)**: Requires HDF5 setup:
```bash
brew install hdf5 c-blosc
export HDF5_DIR=/opt/homebrew/opt/hdf5
export BLOSC_DIR=/opt/homebrew/opt/c-blosc
```

**Windows**: Not officially tested; use WSL2 or Docker

## Data Sources

### BRENDA Database
**Used by**: EnzymeMap

- **Description**: Comprehensive enzyme information system
- **URL**: https://www.brenda-enzymes.org/
- **Size**: 83,000+ enzymes, 13,000+ organisms
- **License**: Free for academic use

### BKMS Database
**Used by**: BKMS Expander

- **Description**: Biochemical pathways and reactions
- **URL**: Previously at biochemical-pathways.com
- **Curated**: Metabolic reactions with organism information

### RetroRules Database
**Used by**: RetroRules Expander

- **Description**: Automated template extraction from metabolic databases
- **Templates**: Multiple diameter levels (2, 4, 6, 8, 10, 12, 14, 16)
- **Biological Score**: Pathway prevalence-based scoring

### Commercial Building Blocks
**Used by**: CommercialSME

- **Description**: Commercially available molecules database
- **Sources**: Aggregated from chemical suppliers
- **Size**: ~1M molecules
- **Updated**: Periodically

### E. coli Metabolism
**Used by**: EcoliSME

- **Description**: E. coli metabolic network
- **Source**: iML1515 model (Monk et al., 2017)
- **Metabolites**: 1,877 metabolites
- **Reactions**: 2,712 reactions

## Example SMILES for Testing

### Simple Molecules
```python
# Alcohols
'CCO'                    # Ethanol
'CCCO'                   # Propanol
'CC(C)CO'                # Isobutanol

# Aldehydes/Ketones
'CCCC=O'                 # Butanal
'CC(C)=O'                # Acetone
'CC(=O)CC'               # Butanone

# Acids
'CC(=O)O'                # Acetic acid
'CCC(=O)O'               # Propionic acid
'CC(O)C(=O)O'            # Lactic acid
```

### Natural Products
```python
# Phenylpropanoids
'c1ccc(cc1)CCO'          # Phenylethanol
'COc1cc(ccc1O)C=O'       # Vanillin
'C=1C=C(C=CC1O)C=CC2=CC(=CC(=C2)O)O'  # Resveratrol

# Flavonoids
'O=C1C(O)=C(Oc2cc(O)cc(O)c12)c1ccc(O)c(O)c1'  # Quercetin

# Alkaloids
'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'  # Caffeine
```

### Pharmaceuticals
```python
# NSAIDs
'CC(=O)Oc1ccccc1C(=O)O'              # Aspirin
'CC(C)Cc1ccc(cc1)C(C)C(=O)O'         # Ibuprofen

# Other drugs
'c1ccc(cc1)C(c2ccccc2)N3CCN(CC3)CCO' # Cetirizine scaffold
```

### Amino Acids
```python
'C(C(=O)O)N'             # Glycine
'CC(C(=O)O)N'            # Alanine
'CC(C)C(C(=O)O)N'        # Valine
'CCC(C)C(C(=O)O)N'       # Isoleucine
```

## Configuration Examples

### Quick Screening Setup
```python
from rbc2 import MCTS, get_expanders

expanders = get_expanders(['retrobiocat'])
mcts = MCTS(target, expanders)
mcts.config.max_search_time = 10
mcts.config.max_depth = 6
```

### Comprehensive Search Setup
```python
from rbc2 import MCTS, get_expanders

expanders = get_expanders([
    'retrobiocat', 'enzymemap', 'bkms',
    'retrorules', 'aizynthfinder'
])

mcts = MCTS(target, expanders)
mcts.config.max_search_time = 300  # 5 minutes
mcts.config.max_depth = 15
mcts.config.max_iterations = 10000
```

### Biocatalysis-Only Setup
```python
from rbc2 import RetroBioCatExpander, EnzymeMapExpander, MCTS

expanders = [
    RetroBioCatExpander(similarity_cutoff=0.6),
    EnzymeMapExpander(cutoff_number=30)
]

mcts = MCTS(target, expanders)
```

## Troubleshooting Guide

### Installation Issues

**Problem**: `ModuleNotFoundError: No module named 'tables'`
```bash
# Solution 1: Install pytables separately
pip install tables

# Solution 2 (ARM Mac): Install with HDF5
brew install hdf5 c-blosc
export HDF5_DIR=/opt/homebrew/opt/hdf5
pip install tables
```

**Problem**: Download failures on first run
```bash
# Ensure internet connection
# Ensure sufficient disk space (5GB)
# Check firewall/proxy settings
# Retry initialization
python -c "from rbc2 import RetroBioCatExpander; RetroBioCatExpander()"
```

### Performance Issues

**Problem**: Very slow searches
```python
# Reduce number of expanders
expanders = get_expanders(['retrobiocat'])  # Instead of all

# Reduce cutoff numbers
expander = EnzymeMapExpander(cutoff_number=20)

# Set strict time limits
mcts.config.max_search_time = 30
```

**Problem**: Memory errors
```python
# Reduce max iterations
mcts.config.max_iterations = 1000

# Process targets in batches
# Use generators instead of loading all at once
```

### No Pathways Found

**Solutions**:
1. Increase search time: `mcts.config.max_search_time = 120`
2. Increase max depth: `mcts.config.max_depth = 15`
3. Use more expanders
4. Lower similarity thresholds
5. Check target molecule validity

## Comparison with Related Tools

| Feature | RetroBioCat 2.0 | READRetro | RetroPath |
|---------|----------------|-----------|-----------|
| **Approach** | Multiple expanders | Dual-view ensemble | Rule-based |
| **Data** | BRENDA, BKMS, custom | KEGG | MetaCyc, KEGG |
| **Chemistry** | Hybrid (bio + chem) | Biosynthesis only | Biosynthesis only |
| **Algorithm** | MCTS | Retro* search | BFS/DFS |
| **Customization** | High (modular) | Medium | Low |
| **Interface** | Python API | Python + Web | Web only |
| **Precedents** | Yes (multiple DBs) | KEGG reactions | Template metadata |

## Best Practices Summary

1. **Start Simple**: Begin with single expander (RetroBioCat)
2. **Expand Gradually**: Add more expanders as needed
3. **Use Appropriate SME**: Commercial for lab work, E. coli for metabolic engineering
4. **Set Realistic Time Limits**: 30-60s for most applications
5. **Validate Results**: Check precedents and availability
6. **Save Pathways**: Export as JSON for reproducibility
7. **Batch Processing**: Use parallel processing for multiple targets

## Additional Resources

### Tutorials and Workshops
- RetroBioCat Documentation: https://retrobiocat-2.readthedocs.io/
- Manchester Biocatalysis Meeting tutorials
- YouTube: Search "RetroBioCat tutorial"

### Related Databases
- **BRENDA**: https://www.brenda-enzymes.org/
- **KEGG**: https://www.genome.jp/kegg/
- **MetaCyc**: https://metacyc.org/
- **Rhea**: https://www.rhea-db.org/

### Related Software
- **RetroPath2.0**: https://github.com/brsynth/RetroPath2
- **RetroPath3.0**: https://github.com/brsynth/rp3
- **Selenzyme**: https://github.com/brsynth/selenzyme
- **EnzymeMap**: https://github.com/hesther/enzymemap

### Literature Reviews
- Finnigan, W., et al. (2020). "Applications of biocatalysis in the pharmaceutical industry"
- Devine, P. N., et al. (2018). "Extending the application of biocatalysis to meet the challenges of drug development"
- Wu, S., et al. (2021). "Machine learning for biosynthetic pathway prediction"

## Citation

If you use RetroBioCat 2.0 in your research, please cite:

```bibtex
@article{finnigan2021retrobiocat,
  title={RetroBioCat as a computer-aided synthesis planning tool for biocatalytic reactions and cascades},
  author={Finnigan, William and Hepworth, Lorna J and Flitsch, Sabine L and Turner, Nicholas J},
  journal={Nature Catalysis},
  volume={4},
  number={1},
  pages={98--104},
  year={2021},
  publisher={Nature Publishing Group}
}
```

For specific expanders, cite the relevant papers listed in the Expander-Specific References section above.

## Version History

- **v2.0 (2023)**: Complete rewrite with modular architecture
  - Multiple expander support
  - MCTS-based multi-step planning
  - Starting material evaluators
  - No web interface (Python package only)

- **v1.0 (2021)**: Original RetroBioCat
  - Web application
  - Manual curation focus
  - Single expansion method

## Support and Community

- **GitHub Issues**: https://github.com/willfinnigan/RetroBioCat-2/issues
- **Documentation**: https://retrobiocat-2.readthedocs.io/
- **Contact**: Check GitHub repository for maintainer contact

## License

RetroBioCat 2.0 is released under the MIT License. See repository for details.

Individual expanders and data sources may have their own licenses - check respective documentation.
