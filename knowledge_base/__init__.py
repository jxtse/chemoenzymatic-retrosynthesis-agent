"""
Unified Chemical-Enzyme-Kinetics Knowledge Base

This package integrates multiple biochemical databases into a single unified API:
- BKMS (Biochemical Kinetic Model System)
- BRENDA (Enzyme kinetics database)
- EnzyExtract (Extracted enzyme data)
- KEGG (Kyoto Encyclopedia of Genes and Genomes)
- USPTO (Patent data)
- UniProt (Protein sequences)
- PubChem (Chemical compounds)
- ZINC (Purchasable compounds)
- RetroBioCat (Biocatalysis retrosynthesis)
"""

from .api import KnowledgeBaseAPI
from .builder import KnowledgeBaseBuilder
from .config import KnowledgeBaseConfig

__all__ = [
    "KnowledgeBaseAPI",
    "KnowledgeBaseBuilder",
    "KnowledgeBaseConfig",
]

__version__ = "0.1.0"
