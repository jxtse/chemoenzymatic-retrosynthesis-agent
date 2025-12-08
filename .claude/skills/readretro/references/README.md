# READRetro Skill References

## Paper

Lee, S., Kim, T., Choi, M.-S., Kwak, Y., Park, J., Hwang, S. J., & Kim, S.-G. (2023). **READRetro: Natural Product Biosynthesis Planning with Retrieval-Augmented Dual-View Retrosynthesis**. *bioRxiv*.

https://www.biorxiv.org/content/10.1101/2023.03.21.533616v1

## Official Resources

- **Web Version**: https://readretro.net
- **GitHub Repository**: https://github.com/SeulLee05/READRetro
- **Data (Zenodo)**: https://zenodo.org/records/11485641

## Key Datasets

1. **KEGG PATHWAY Database**
   - Metabolic pathway database
   - https://www.genome.jp/kegg/pathway.html

2. **KEGG REACTION Database**
   - Enzymatic reaction templates
   - https://www.genome.jp/kegg/reaction/

3. **Biochemical Dataset**
   - Natural product-specific retrosynthesis
   - Included in Zenodo data

## Related Models

### Retroformer
- Template-based retrosynthesis model
- Paper: Wan, Y., et al. (2022). "Retroformer: Pushing the Limits of Interpretation in End-to-End Retrosynthesis with Transformer and Graph Transformer"
- GitHub: https://github.com/yuewan2/Retroformer

### Graph2SMILES
- Template-free graph-to-sequence model
- Paper: Tu, Z., & Coley, C. W. (2022). "Permutation invariant graph-to-sequence model for template-free retrosynthesis and reaction prediction"
- GitHub: https://github.com/coleygroup/Graph2SMILES

## Installation Dependencies

### Core Requirements
- Python 3.8
- PyTorch 1.12.0 (CUDA 11.3 compatible)
- RDKit
- OpenNMT-py 2.3.0

### Full Dependency List
```
torch==1.12.0+cu113
torchvision==0.13.0+cu113
rdkit-pypi
OpenNMT-py==2.3.0
networkx==2.5
numpy==1.22
pandas
scipy
gin-config==0.4.0
easydict
tqdm
```

## System Requirements

### Hardware
- **Minimum**: 8GB RAM, CPU only (slow)
- **Recommended**: 16GB RAM, NVIDIA GPU with 8GB+ VRAM
- **Optimal**: 32GB RAM, NVIDIA GPU with 12GB+ VRAM

### Storage
- Model files: ~600 MB
- KEGG data: ~20 MB
- Working space: 1 GB recommended

### Operating System
- Linux (tested on CentOS 7)
- Ubuntu 18.04+
- macOS (CPU only)

## Environment Location

**Conda Environment**: `/home/xiejinxiang/.conda/envs/readretro`

**READRetro Directory**: `/home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/READRetro`

**Skill Directory**: `/home/xiejinxiang/chemoenzymatic-retrosynthesis-agent/.claude/skills/readretro`

## Common SMILES Examples

### Natural Products

```python
# Resveratrol (stilbenoid)
'C=1C=C(C=CC1O)C=CC2=CC(=CC(=C2)O)O'

# Curcumin (diarylheptanoid)
'COc1cc(ccc1O)C=CC(=O)CC(=O)C=Cc1ccc(O)c(OC)c1'

# Quercetin (flavonoid)
'O=C1C(O)=C(Oc2cc(O)cc(O)c12)c1ccc(O)c(O)c1'

# Caffeine (alkaloid)
'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
```

### Terpenes

```python
# Geraniol (monoterpene)
'CC(C)=CCCC(C)=CCO'

# Limonene (monoterpene)
'CC(=C)C1CCC(=C)CC1'

# Farnesol (sesquiterpene)
'CC(C)=CCCC(C)=CCCC(C)=CCO'
```

### Polyketides

```python
# 6-Methylsalicylic acid
'Cc1cccc(O)c1C(=O)O'

# Orsellinic acid
'Cc1cc(O)cc(O)c1C(=O)O'
```

## Troubleshooting Reference

### Error Messages

1. **"ModuleNotFoundError: No module named 'scipy'"**
   - Solution: `pip install scipy`

2. **"ModuleNotFoundError: No module named 'gin'"**
   - Solution: `pip install gin-config==0.4.0`

3. **"CUDA out of memory"**
   - Solutions:
     - Reduce beam_size
     - Reduce num_threads in run_mp.py
     - Use CPU mode (slower)

4. **"Invalid SMILES"**
   - Solution: Validate with RDKit first
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles(smiles)
   canonical = Chem.MolToSmiles(mol)
   ```

5. **"configurable() got an unexpected keyword argument 'whitelist'"**
   - Solution: Downgrade gin-config
   ```bash
   pip install gin-config==0.4.0
   ```

### Performance Tips

1. **GPU Selection**: Use newest GPU available
   ```bash
   CUDA_VISIBLE_DEVICES=0 python run.py 'SMILES'
   ```

2. **Memory Optimization**: Reduce batch/beam size
   ```python
   # In run_mp.py
   num_threads = 2  # Instead of 4
   beam_size = 5    # Instead of 10
   ```

3. **Speed vs Accuracy**: Adjust search parameters
   - Fast: beam_size=5, num_path=50
   - Balanced: beam_size=10, num_path=100
   - Thorough: beam_size=20, num_path=200

## Citation

If you use READRetro in your research, please cite:

```bibtex
@article{lee2023readretro,
  title={READRetro: Natural Product Biosynthesis Planning with Retrieval-Augmented Dual-View Retrosynthesis},
  author={Lee, Seul and Kim, Taein and Choi, Min-Soo and Kwak, Yejin and Park, Jeongbin and Hwang, Sung Ju and Kim, Sang-Gyu},
  journal={bioRxiv},
  year={2023},
  url={https://www.biorxiv.org/content/10.1101/2023.03.21.533616v1}
}
```

## Additional Resources

### Natural Product Databases
- **Natural Products Atlas**: https://www.npatlas.org/
- **PubChem**: https://pubchem.ncbi.nlm.nih.gov/
- **ChEMBL**: https://www.ebi.ac.uk/chembl/

### Biosynthesis Tools
- **antiSMASH**: Secondary metabolite biosynthesis analysis
- **RetroPath2.0**: Metabolic pathway design
- **BNICE.ch**: Biochemical network synthesis

### Learning Resources
- **KEGG Pathway Tutorial**: https://www.genome.jp/kegg/pathway.html
- **Retrosynthesis Review**: Coley, C. W., et al. (2018). "A graph-convolutional neural network model for the prediction of chemical reactivity"
- **Natural Product Biosynthesis**: Dewick, P. M. (2009). "Medicinal Natural Products: A Biosynthetic Approach"
