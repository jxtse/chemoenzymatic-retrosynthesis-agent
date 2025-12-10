"""
Utility functions for the retrosynthesis agent
"""

import json
import logging
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


def name_to_smiles(compound_name: str) -> str:
    """
    Convert compound name to SMILES string.

    使用PubChem查询分子名称对应的SMILES

    Args:
        compound_name: 化合物名称 (中英文均可)

    Returns:
        JSON string with SMILES and basic info

    Example:
        result = name_to_smiles("布洛芬")
        result = name_to_smiles("ibuprofen")
    """
    try:
        import pubchempy as pcp

        # 查询PubChem
        compounds = pcp.get_compounds(compound_name, 'name')

        if not compounds:
            return json.dumps({
                "success": False,
                "query": compound_name,
                "error": f"没找到 '{compound_name}'，请检查拼写或试试英文名"
            })

        compound = compounds[0]

        result = {
            "success": True,
            "query": compound_name,
            "smiles": compound.isomeric_smiles or compound.canonical_smiles,
            "iupac_name": compound.iupac_name,
            "molecular_formula": compound.molecular_formula,
            "molecular_weight": compound.molecular_weight,
            "cid": compound.cid,
        }

        return json.dumps(result, indent=2, ensure_ascii=False)

    except ImportError:
        return json.dumps({
            "success": False,
            "error": "需要安装pubchempy: uv add pubchempy"
        })
    except Exception as e:
        logger.error(f"Error converting name to SMILES: {e}")
        return json.dumps({
            "success": False,
            "query": compound_name,
            "error": str(e)
        })


def analyze_molecule_properties(smiles: str) -> str:
    """
    Analyze molecular properties using RDKit.

    计算分子的基本性质

    Args:
        smiles: SMILES字符串

    Returns:
        JSON string with molecular properties

    Example:
        result = analyze_molecule_properties("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            return json.dumps({
                "success": False,
                "error": f"无效的SMILES: {smiles}"
            })

        properties = {
            "success": True,
            "smiles": smiles,
            "properties": {
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "logp": round(Descriptors.MolLogP(mol), 2),
                "tpsa": round(Descriptors.TPSA(mol), 2),
                "num_h_donors": Lipinski.NumHDonors(mol),
                "num_h_acceptors": Lipinski.NumHAcceptors(mol),
                "num_rotatable_bonds": Lipinski.NumRotatableBonds(mol),
                "num_aromatic_rings": Lipinski.NumAromaticRings(mol),
                "num_atoms": mol.GetNumAtoms(),
                "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            },
            "drug_likeness": {
                "lipinski_rule_of_5": _check_lipinski(mol),
                "complexity": _estimate_complexity(mol),
            }
        }

        return json.dumps(properties, indent=2, ensure_ascii=False)

    except ImportError:
        return json.dumps({
            "success": False,
            "error": "需要安装rdkit: uv add rdkit"
        })
    except Exception as e:
        logger.error(f"Error analyzing molecule: {e}")
        return json.dumps({
            "success": False,
            "error": str(e)
        })


def _check_lipinski(mol) -> Dict[str, Any]:
    """
    Check Lipinski's Rule of Five.

    Args:
        mol: RDKit mol object

    Returns:
        Dictionary with rule compliance
    """
    from rdkit.Chem import Descriptors, Lipinski

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = 0
    details = {}

    if mw > 500:
        violations += 1
        details["molecular_weight"] = f"❌ {mw:.1f} > 500"
    else:
        details["molecular_weight"] = f"✅ {mw:.1f} ≤ 500"

    if logp > 5:
        violations += 1
        details["logp"] = f"❌ {logp:.1f} > 5"
    else:
        details["logp"] = f"✅ {logp:.1f} ≤ 5"

    if hbd > 5:
        violations += 1
        details["h_donors"] = f"❌ {hbd} > 5"
    else:
        details["h_donors"] = f"✅ {hbd} ≤ 5"

    if hba > 10:
        violations += 1
        details["h_acceptors"] = f"❌ {hba} > 10"
    else:
        details["h_acceptors"] = f"✅ {hba} ≤ 10"

    return {
        "passes": violations <= 1,
        "violations": violations,
        "details": details,
        "summary": "类药性好 ✅" if violations <= 1 else f"类药性差 ({violations}条违反) ❌"
    }


def _estimate_complexity(mol) -> str:
    """
    Estimate synthetic complexity.

    Args:
        mol: RDKit mol object

    Returns:
        Complexity assessment string
    """
    from rdkit import Chem
    from rdkit.Chem import Lipinski

    # 简单的启发式评估
    num_rings = Lipinski.NumAromaticRings(mol) + Lipinski.NumAliphaticRings(mol)
    num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_rot_bonds = Lipinski.NumRotatableBonds(mol)
    num_heavy = mol.GetNumHeavyAtoms()

    score = 0

    # 评分
    if num_heavy < 15:
        score += 1
        complexity = "简单"
    elif num_heavy < 30:
        score += 2
        complexity = "中等"
    else:
        score += 3
        complexity = "复杂"

    if num_rings > 2:
        score += 1

    if num_stereo > 2:
        score += 1

    if score <= 2:
        return f"{complexity} (合成难度低 ⭐)"
    elif score <= 4:
        return f"{complexity} (合成难度中 ⭐⭐)"
    else:
        return f"{complexity} (合成难度高 ⭐⭐⭐)"


def smiles_to_inchi(smiles: str) -> str:
    """
    Convert SMILES to InChI.

    Args:
        smiles: SMILES string

    Returns:
        JSON with InChI and InChIKey
    """
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return json.dumps({"success": False, "error": "Invalid SMILES"})

        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)

        return json.dumps({
            "success": True,
            "smiles": smiles,
            "inchi": inchi,
            "inchikey": inchikey
        }, indent=2)

    except Exception as e:
        return json.dumps({"success": False, "error": str(e)})


def calculate_feasibility_score(
    pathway_data: Dict[str, Any],
    availability_data: Dict[str, Any],
    enzyme_data: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Calculate overall pathway feasibility score.

    综合计算路径可行性评分 (0-10分)

    Args:
        pathway_data: RetroBioCat pathway result
        availability_data: Commercial availability result
        enzyme_data: Optional enzyme information

    Returns:
        Detailed scoring breakdown
    """
    score = 0
    max_score = 10
    breakdown = {}

    # 1. 步数评分 (0-3分)
    num_steps = pathway_data.get("total_steps", 999)
    if num_steps <= 3:
        step_score = 3
        breakdown["steps"] = "✅ 3步以内，很简洁"
    elif num_steps <= 6:
        step_score = 2
        breakdown["steps"] = "⚠️ 4-6步，中等复杂度"
    else:
        step_score = 1
        breakdown["steps"] = "❌ 7步以上，较复杂"
    score += step_score

    # 2. 原料可获得性 (0-3分)
    total_materials = availability_data.get("total_molecules", 0)
    available_count = availability_data.get("available_count", 0)

    if total_materials > 0:
        availability_ratio = available_count / total_materials
        if availability_ratio >= 1.0:
            avail_score = 3
            breakdown["availability"] = "✅ 所有起始原料都能买到"
        elif availability_ratio >= 0.5:
            avail_score = 2
            breakdown["availability"] = f"⚠️ {available_count}/{total_materials} 原料可买到"
        else:
            avail_score = 1
            breakdown["availability"] = f"❌ 仅 {available_count}/{total_materials} 原料可买到"
    else:
        avail_score = 1
        breakdown["availability"] = "⚠️ 未检查原料可获得性"
    score += avail_score

    # 3. 酶的可获得性 (0-2分)
    bio_steps = pathway_data.get("biocatalysis_steps", 0)
    if bio_steps > 0:
        if enzyme_data and enzyme_data.get("found", 0) > 0:
            enzyme_score = 2
            breakdown["enzymes"] = f"✅ 找到了 {bio_steps} 个生物催化步骤的酶"
        else:
            enzyme_score = 1
            breakdown["enzymes"] = f"⚠️ {bio_steps} 个生物催化步骤，需要确认酶可获得性"
    else:
        enzyme_score = 2
        breakdown["enzymes"] = "纯化学路线，不需要酶"
    score += enzyme_score

    # 4. 文献先例 (0-2分)
    reactions = pathway_data.get("reactions", [])
    reactions_with_precedents = sum(
        1 for r in reactions if r.get("precedents", 0) > 0
    )
    if len(reactions) > 0:
        precedent_ratio = reactions_with_precedents / len(reactions)
        if precedent_ratio >= 0.8:
            precedent_score = 2
            breakdown["precedents"] = f"✅ {reactions_with_precedents}/{len(reactions)} 步有文献先例"
        elif precedent_ratio >= 0.5:
            precedent_score = 1
            breakdown["precedents"] = f"⚠️ {reactions_with_precedents}/{len(reactions)} 步有文献先例"
        else:
            precedent_score = 0
            breakdown["precedents"] = f"❌ 仅 {reactions_with_precedents}/{len(reactions)} 步有文献先例"
    else:
        precedent_score = 0
        breakdown["precedents"] = "⚠️ 未找到反应信息"
    score += precedent_score

    # 总分和建议
    if score >= 8:
        recommendation = "💚 强烈推荐！这条路线很靠谱，可以尝试"
        confidence = "高"
    elif score >= 5:
        recommendation = "💛 可以尝试，但需要注意一些难点"
        confidence = "中"
    else:
        recommendation = "❤️ 不太推荐，建议重新设计或换个思路"
        confidence = "低"

    return {
        "total_score": score,
        "max_score": max_score,
        "percentage": round(score / max_score * 100, 1),
        "breakdown": breakdown,
        "recommendation": recommendation,
        "confidence": confidence,
    }
