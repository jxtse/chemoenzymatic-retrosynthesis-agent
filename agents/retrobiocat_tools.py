"""
RetroBioCat Integration Tools for Chemoenzymatic Agent

Provides retrosynthesis planning using RetroBioCat 2.0 framework.
Combines biocatalysis and traditional chemistry approaches.
"""

from typing import List, Dict, Optional, Tuple, Any
import json


class RetroBioCatTools:
    """
    RetroBioCat 2.0 tools for retrosynthesis planning.

    Provides:
    - Single-step retrosynthesis with multiple expanders
    - Multi-step pathway planning with MCTS
    - Starting material availability checking
    - Hybrid bio + chem route design
    """

    def __init__(self, default_expanders: Optional[List[str]] = None):
        """
        Initialize RetroBioCat tools.

        Args:
            default_expanders: Default expanders to use.
                             If None, uses ['retrobiocat', 'enzymemap', 'aizynthfinder']
        """
        self._rbc2_available = None
        self._check_availability()

        # 使用已验证可用的 expanders (包括 BKMS)
        if default_expanders is None:
            self.default_expanders = ['retrobiocat', 'enzymemap', 'bkms', 'aizynthfinder']
        else:
            self.default_expanders = default_expanders

    def _check_availability(self):
        """Check if RetroBioCat 2.0 is installed."""
        try:
            import rbc2
            self._rbc2_available = True
        except ImportError:
            self._rbc2_available = False

    def _ensure_installed(self):
        """Ensure RetroBioCat is installed."""
        if not self._rbc2_available:
            raise ImportError(
                "RetroBioCat 2.0 is not installed. "
                "Install with: uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git"
            )

    def plan_biocatalytic_route(
        self,
        target_smiles: str,
        max_steps: int = 5,
        use_chemistry: bool = True,
        max_search_time: int = 30,
        custom_expanders: Optional[List[str]] = None
    ) -> str:
        """
        Plan biocatalytic synthesis route using MCTS.

        Args:
            target_smiles: Target molecule SMILES string
            max_steps: Maximum number of synthetic steps
            use_chemistry: Include traditional chemistry steps (hybrid routes)
            max_search_time: Maximum search time in seconds
            custom_expanders: Custom list of expanders to use (overrides use_chemistry)

        Returns:
            JSON string with pathway results
        """
        self._ensure_installed()

        from rbc2 import MCTS, get_expanders, CommercialSME

        # Select expanders
        if custom_expanders:
            expander_names = custom_expanders
        elif use_chemistry:
            # 混合路线：生物催化 + 化学反应（使用已验证的 expanders，包含 BKMS）
            expander_names = ['retrobiocat', 'enzymemap', 'bkms', 'aizynthfinder']
        else:
            # 纯生物催化（包含 BKMS）
            expander_names = ['retrobiocat', 'enzymemap', 'bkms']

        expanders = get_expanders(expander_names)
        sme = CommercialSME()

        # Run MCTS
        mcts = MCTS(target_smiles, expanders, starting_material_evaluator=sme)
        mcts.config.max_search_time = max_search_time
        mcts.config.max_depth = max_steps
        mcts.run()

        # Get results
        solved_pathways = mcts.get_solved_pathways()

        if not solved_pathways:
            return json.dumps({
                "target": target_smiles,
                "pathways_found": 0,
                "message": "No complete pathways found. Try increasing max_search_time or max_steps."
            })

        # Format results
        results = {
            "target": target_smiles,
            "pathways_found": len(solved_pathways),
            "pathways": []
        }

        for i, pathway in enumerate(solved_pathways[:5]):  # Top 5
            bio_steps = sum(1 for r in pathway.reactions if r.rxn_domain == 'biocatalysis')
            chem_steps = sum(1 for r in pathway.reactions if r.rxn_domain == 'chemistry')

            pathway_data = {
                "pathway_id": i + 1,
                "total_steps": pathway.pathway_length,
                "biocatalysis_steps": bio_steps,
                "chemistry_steps": chem_steps,
                "reactions": [],
                "starting_materials": []
            }

            # Reaction sequence
            for j, rxn in enumerate(pathway.reactions):
                reaction_data = {
                    "step": j + 1,
                    "reaction_smiles": rxn.reaction_smiles(),
                    "name": rxn.name,
                    "type": rxn.rxn_type,
                    "domain": rxn.rxn_domain,
                    "score": round(rxn.score, 3)
                }

                if rxn.precedents:
                    reaction_data["precedents"] = len(rxn.precedents)
                    reaction_data["top_precedent"] = {
                        "name": rxn.precedents[0].name,
                        "similarity": round(rxn.precedents[0].similarity, 3)
                    }

                pathway_data["reactions"].append(reaction_data)

            # Starting materials
            for smi in pathway.end_smis():
                available, info = sme.eval(smi)
                pathway_data["starting_materials"].append({
                    "smiles": smi,
                    "commercially_available": available,
                    "info": info
                })

            results["pathways"].append(pathway_data)

        return json.dumps(results, indent=2)

    def find_enzymatic_reactions(
        self,
        target_smiles: str,
        expander_types: Optional[List[str]] = None
    ) -> str:
        """
        Find one-step enzymatic reactions for target molecule.

        Args:
            target_smiles: Target molecule SMILES
            expander_types: List of expander types to use
                          Options: 'retrobiocat', 'enzymemap', 'aizynthfinder', 'bkms', 'retrorules'
                          Default: ['retrobiocat', 'enzymemap'] (bio-only)

        Returns:
            JSON string with reaction results
        """
        self._ensure_installed()

        from rbc2 import get_expanders

        if expander_types is None:
            # 默认只用生物催化（不包含化学反应，包含 BKMS）
            expander_types = ['retrobiocat', 'enzymemap', 'bkms']

        expanders = get_expanders(expander_types)

        all_reactions = []
        expander_results = {}

        for expander in expanders:
            reactions = expander.get_reactions(target_smiles)
            expander_name = expander.__class__.__name__
            expander_results[expander_name] = len(reactions)
            all_reactions.extend(reactions)

        # Format results
        results = {
            "target": target_smiles,
            "total_reactions": len(all_reactions),
            "by_expander": expander_results,
            "reactions": []
        }

        # Sort by score
        all_reactions.sort(key=lambda r: r.score, reverse=True)

        for i, rxn in enumerate(all_reactions[:20]):  # Top 20
            reaction_data = {
                "rank": i + 1,
                "reaction_smiles": rxn.reaction_smiles(),
                "name": rxn.name,
                "score": round(rxn.score, 3),
                "type": rxn.rxn_type,
                "domain": rxn.rxn_domain,
                "substrates": rxn.substrates
            }

            # Literature precedents
            if rxn.precedents:
                precedent_data = []
                for prec in rxn.precedents[:3]:  # Top 3 precedents
                    precedent_data.append({
                        "name": prec.name,
                        "similarity": round(prec.similarity, 3),
                        "data": prec.data
                    })
                reaction_data["precedents"] = precedent_data

            # Template info
            if hasattr(rxn, 'template_metadata') and rxn.template_metadata:
                reaction_data["template"] = rxn.template_metadata

            results["reactions"].append(reaction_data)

        return json.dumps(results, indent=2)

    def check_commercial_availability(
        self,
        smiles_list: List[str]
    ) -> str:
        """
        Check commercial availability of molecules.

        Args:
            smiles_list: List of SMILES strings to check

        Returns:
            JSON string with availability results
        """
        self._ensure_installed()

        from rbc2 import CommercialSME

        sme = CommercialSME()

        results = {
            "total_molecules": len(smiles_list),
            "available_count": 0,
            "molecules": []
        }

        for smi in smiles_list:
            available, info = sme.eval(smi)

            if available:
                results["available_count"] += 1

            results["molecules"].append({
                "smiles": smi,
                "available": available,
                "info": info
            })

        return json.dumps(results, indent=2)

    def compare_retrosynthesis_approaches(
        self,
        target_smiles: str
    ) -> str:
        """
        Compare different retrosynthesis approaches for target.

        Compares:
        - RetroBioCat (curated enzymatic)
        - EnzymeMap (BRENDA-based)
        - BKMS (metabolic)
        - AIZynthfinder (chemistry)

        Args:
            target_smiles: Target molecule SMILES

        Returns:
            JSON string comparing approaches
        """
        self._ensure_installed()

        from rbc2 import (RetroBioCatExpander, EnzymeMapExpander,
                         BKMSExpander, AIZynthfinderExpander)

        approaches = {
            "RetroBioCat": RetroBioCatExpander(),
            "EnzymeMap": EnzymeMapExpander(),
            "BKMS": BKMSExpander(),
            "AIZynthfinder": AIZynthfinderExpander()
        }

        results = {
            "target": target_smiles,
            "approaches": {}
        }

        for name, expander in approaches.items():
            try:
                reactions = expander.get_reactions(target_smiles)

                approach_data = {
                    "reactions_found": len(reactions),
                    "top_reactions": []
                }

                for rxn in reactions[:5]:  # Top 5
                    approach_data["top_reactions"].append({
                        "reaction": rxn.reaction_smiles(),
                        "score": round(rxn.score, 3),
                        "domain": rxn.rxn_domain,
                        "has_precedents": len(rxn.precedents) > 0 if rxn.precedents else False
                    })

                results["approaches"][name] = approach_data

            except Exception as e:
                results["approaches"][name] = {
                    "error": str(e)
                }

        return json.dumps(results, indent=2)

    def analyze_pathway_feasibility(
        self,
        target_smiles: str,
        pathway_id: int = 1,
        max_search_time: int = 30
    ) -> str:
        """
        Analyze feasibility of best retrosynthesis pathway.

        Args:
            target_smiles: Target molecule
            pathway_id: Which pathway to analyze (1-based)
            max_search_time: Search time limit

        Returns:
            JSON string with detailed pathway analysis
        """
        self._ensure_installed()

        from rbc2 import MCTS, get_expanders, CommercialSME
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        expanders = get_expanders(['retrobiocat', 'enzymemap', 'bkms', 'aizynthfinder'])
        sme = CommercialSME()

        mcts = MCTS(target_smiles, expanders, starting_material_evaluator=sme)
        mcts.config.max_search_time = max_search_time
        mcts.run()

        pathways = mcts.get_solved_pathways()

        if not pathways:
            return json.dumps({"error": "No pathways found"})

        if pathway_id > len(pathways):
            pathway_id = 1

        pathway = pathways[pathway_id - 1]

        # Detailed analysis
        analysis = {
            "pathway_id": pathway_id,
            "target": {
                "smiles": target_smiles,
                "molecular_weight": 0,
                "logp": 0
            },
            "pathway_metrics": {
                "total_steps": pathway.pathway_length,
                "total_reactions": len(pathway.reactions),
                "biocatalysis_steps": 0,
                "chemistry_steps": 0
            },
            "starting_materials": {
                "count": 0,
                "all_available": True,
                "details": []
            },
            "reaction_details": [],
            "feasibility_score": 0.0
        }

        # Target properties
        try:
            mol = Chem.MolFromSmiles(target_smiles)
            analysis["target"]["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
            analysis["target"]["logp"] = round(Descriptors.MolLogP(mol), 2)
        except:
            pass

        # Count domains
        for rxn in pathway.reactions:
            if rxn.rxn_domain == 'biocatalysis':
                analysis["pathway_metrics"]["biocatalysis_steps"] += 1
            elif rxn.rxn_domain == 'chemistry':
                analysis["pathway_metrics"]["chemistry_steps"] += 1

        # Starting materials analysis
        end_smis = pathway.end_smis()
        analysis["starting_materials"]["count"] = len(end_smis)

        for smi in end_smis:
            available, info = sme.eval(smi)
            analysis["starting_materials"]["details"].append({
                "smiles": smi,
                "available": available,
                "source": info.get("source", "unknown") if info else "unknown"
            })

            if not available:
                analysis["starting_materials"]["all_available"] = False

        # Reaction details with precedents
        for i, rxn in enumerate(pathway.reactions):
            rxn_detail = {
                "step": i + 1,
                "reaction": rxn.reaction_smiles(),
                "name": rxn.name,
                "domain": rxn.rxn_domain,
                "score": round(rxn.score, 3),
                "has_precedents": False,
                "precedent_count": 0
            }

            if rxn.precedents:
                rxn_detail["has_precedents"] = True
                rxn_detail["precedent_count"] = len(rxn.precedents)
                rxn_detail["best_precedent_similarity"] = round(rxn.precedents[0].similarity, 3)

            analysis["reaction_details"].append(rxn_detail)

        # Calculate feasibility score (0-1)
        feasibility = 1.0

        # Penalize long pathways
        if pathway.pathway_length > 5:
            feasibility *= 0.8

        # Penalize unavailable starting materials
        if not analysis["starting_materials"]["all_available"]:
            feasibility *= 0.6

        # Reward precedents
        precedent_ratio = sum(1 for r in pathway.reactions if r.precedents) / len(pathway.reactions)
        feasibility *= (0.5 + 0.5 * precedent_ratio)

        analysis["feasibility_score"] = round(feasibility, 3)

        return json.dumps(analysis, indent=2)
