"""
Knowledge Base Tools for Autogen Agents

Provides tool functions that wrap the unified knowledge base API
for use by Autogen agents in retrieval-augmented generation.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from knowledge_base import KnowledgeBaseAPI

logger = logging.getLogger(__name__)


class KnowledgeBaseTools:
    """
    Wrapper class providing tool functions for Autogen agents to access
    the unified chemical-enzyme-kinetics knowledge base.
    """

    def __init__(self, kb_path: str | Path):
        """
        Initialize knowledge base tools.

        Args:
            kb_path: Path to the knowledge base file (JSONL or Parquet)
        """
        self.kb_path = Path(kb_path)
        self.api = KnowledgeBaseAPI(self.kb_path)
        self._loaded = False

    def _ensure_loaded(self) -> None:
        """Ensure knowledge base is loaded."""
        if not self._loaded:
            logger.info(f"Loading knowledge base from {self.kb_path}...")
            self.api.load()
            self._loaded = True
            logger.info("Knowledge base loaded successfully")

    def search_enzyme_by_ec(
        self,
        ec_number: str,
        include_sources: Optional[List[str]] = None,
    ) -> str:
        """
        Search for enzyme information by EC number.

        Args:
            ec_number: EC number (e.g., "1.1.1.1")
            include_sources: Optional list of data sources to filter
                           (e.g., ["BRENDA", "BKMS"])

        Returns:
            JSON string with search results including:
            - Enzyme names
            - Organisms
            - Reaction equations
            - Kinetic parameters (kcat, Km)
            - Sequences (if available)

        Example:
            result = search_enzyme_by_ec("1.1.1.1")
        """
        self._ensure_loaded()

        try:
            result = self.api.query_by_ec(ec_number, include_sources=include_sources)

            # Format for LLM consumption
            formatted = {
                "query": {"ec_number": ec_number},
                "found": result["count"],
                "enzymes": [],
            }

            for record in result["records"][:10]:  # Limit to 10 for brevity
                enzyme_info = {
                    "id": record.get("id"),
                    "source": record.get("source", {}).get("dataset"),
                    "enzyme_name": record.get("enzyme", {}).get("name"),
                    "organism": record.get("enzyme", {}).get("organism"),
                    "reaction": record.get("reaction", {}).get("equation_text"),
                    "kcat": record.get("kinetics", {}).get("kcat", {}).get("value"),
                    "km": record.get("kinetics", {}).get("km", {}).get("value"),
                    "has_sequence": bool(record.get("enzyme", {}).get("sequence")),
                    "uniprot_ids": record.get("enzyme", {}).get("uniprot_ids", []),
                }
                formatted["enzymes"].append(enzyme_info)

            return json.dumps(formatted, indent=2)

        except Exception as e:
            logger.error(f"Error searching by EC: {e}")
            return json.dumps({"error": str(e), "query": ec_number})

    def search_reactions_by_compound(
        self,
        compound_name: str,
        role: Optional[str] = None,
    ) -> str:
        """
        Search for reactions involving a specific compound.

        Args:
            compound_name: Name of the compound (e.g., "glucose")
            role: Optional role filter - "substrate" or "product"

        Returns:
            JSON string with reactions involving the compound

        Example:
            result = search_reactions_by_compound("glucose", role="substrate")
        """
        self._ensure_loaded()

        try:
            result = self.api.query_by_compound(compound_name, role=role)

            formatted = {
                "query": {"compound": compound_name, "role": role},
                "found": result["count"],
                "reactions": [],
            }

            for record in result["records"][:10]:
                reaction_info = {
                    "id": record.get("id"),
                    "ec_number": record.get("primary_ec"),
                    "enzyme_name": record.get("enzyme", {}).get("name"),
                    "reaction": record.get("reaction", {}).get("equation_text"),
                    "direction": record.get("reaction", {}).get("direction"),
                    "substrates": [
                        s.get("name")
                        for s in record.get("reaction", {}).get("substrates", [])
                    ],
                    "products": [
                        p.get("name")
                        for p in record.get("reaction", {}).get("products", [])
                    ],
                    "organism": record.get("enzyme", {}).get("organism"),
                }
                formatted["reactions"].append(reaction_info)

            return json.dumps(formatted, indent=2)

        except Exception as e:
            logger.error(f"Error searching by compound: {e}")
            return json.dumps({"error": str(e), "query": compound_name})

    def find_retrosynthesis_pathway(
        self,
        target_compound: str,
        available_substrates: Optional[List[str]] = None,
    ) -> str:
        """
        Find possible retrosynthesis pathways for a target compound.

        Args:
            target_compound: Target compound to synthesize
            available_substrates: List of available starting materials

        Returns:
            JSON string with possible enzymatic pathways

        Example:
            result = find_retrosynthesis_pathway(
                "glucose-6-phosphate",
                available_substrates=["glucose", "ATP"]
            )
        """
        self._ensure_loaded()

        try:
            # Find reactions that produce the target
            products_result = self.api.query_by_compound(
                target_compound, role="product"
            )

            pathways = []

            for record in products_result["records"][:10]:
                reaction = record.get("reaction", {})
                substrates = [s.get("name") for s in reaction.get("substrates", [])]

                # Check if substrates are available
                pathway_feasible = True
                if available_substrates:
                    pathway_feasible = all(
                        any(sub.lower() in avail.lower() for avail in available_substrates)
                        for sub in substrates
                    )

                pathway = {
                    "ec_number": record.get("primary_ec"),
                    "enzyme": record.get("enzyme", {}).get("name"),
                    "reaction": reaction.get("equation_text"),
                    "substrates_required": substrates,
                    "feasible": pathway_feasible,
                    "kcat": record.get("kinetics", {}).get("kcat", {}).get("value"),
                    "km": record.get("kinetics", {}).get("km", {}).get("value"),
                    "organism": record.get("enzyme", {}).get("organism"),
                }
                pathways.append(pathway)

            result = {
                "target": target_compound,
                "available_substrates": available_substrates or [],
                "pathways_found": len(pathways),
                "pathways": pathways,
            }

            return json.dumps(result, indent=2)

        except Exception as e:
            logger.error(f"Error finding pathway: {e}")
            return json.dumps({"error": str(e), "target": target_compound})

    def get_kinetic_parameters(
        self,
        ec_number: Optional[str] = None,
        substrate: Optional[str] = None,
        organism: Optional[str] = None,
    ) -> str:
        """
        Retrieve kinetic parameters (kcat, Km) with optional filters.

        Args:
            ec_number: Filter by EC number
            substrate: Filter by substrate name
            organism: Filter by organism name

        Returns:
            JSON string with kinetic parameters

        Example:
            result = get_kinetic_parameters(
                ec_number="1.1.1.1",
                organism="Escherichia coli"
            )
        """
        self._ensure_loaded()

        try:
            result = self.api.get_kinetics(
                ec_number=ec_number,
                substrate=substrate,
                organism=organism,
            )

            formatted = {
                "query": {
                    "ec_number": ec_number,
                    "substrate": substrate,
                    "organism": organism,
                },
                "measurements_found": result["count"],
                "kinetics": result["kinetics"][:20],  # Limit to 20
            }

            return json.dumps(formatted, indent=2)

        except Exception as e:
            logger.error(f"Error getting kinetics: {e}")
            return json.dumps({"error": str(e)})

    def get_enzyme_sequence(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> str:
        """
        Retrieve enzyme amino acid sequences by EC number.

        Args:
            ec_number: EC number
            organism: Optional organism filter

        Returns:
            JSON string with enzyme sequences and UniProt IDs

        Example:
            result = get_enzyme_sequence("1.1.1.1", organism="E. coli")
        """
        self._ensure_loaded()

        try:
            result = self.api.query_by_ec(ec_number)

            sequences = []
            for record in result["records"]:
                enzyme = record.get("enzyme", {})
                org = enzyme.get("organism", "")

                # Filter by organism if specified
                if organism and organism.lower() not in org.lower():
                    continue

                sequence = enzyme.get("sequence")
                if sequence:
                    sequences.append({
                        "organism": org,
                        "sequence": sequence[:100] + "..." if len(sequence) > 100 else sequence,
                        "full_length": len(sequence),
                        "uniprot_ids": enzyme.get("uniprot_ids", []),
                        "pdb_ids": enzyme.get("pdb_ids", []),
                    })

            formatted = {
                "ec_number": ec_number,
                "organism_filter": organism,
                "sequences_found": len(sequences),
                "sequences": sequences[:5],  # Limit to 5
            }

            return json.dumps(formatted, indent=2)

        except Exception as e:
            logger.error(f"Error getting sequences: {e}")
            return json.dumps({"error": str(e)})

    def compare_enzymes(
        self,
        ec_numbers: List[str],
    ) -> str:
        """
        Compare multiple enzymes by their kinetic parameters.

        Args:
            ec_numbers: List of EC numbers to compare

        Returns:
            JSON string with comparison data

        Example:
            result = compare_enzymes(["1.1.1.1", "1.1.1.2"])
        """
        self._ensure_loaded()

        try:
            comparison = []

            for ec in ec_numbers:
                result = self.api.query_by_ec(ec)

                if result["records"]:
                    # Aggregate kinetics
                    kcat_values = []
                    km_values = []

                    for record in result["records"]:
                        kinetics = record.get("kinetics", {})
                        kcat = kinetics.get("kcat", {}).get("value")
                        km = kinetics.get("km", {}).get("value")

                        if kcat:
                            kcat_values.append(kcat)
                        if km:
                            km_values.append(km)

                    comparison.append({
                        "ec_number": ec,
                        "enzyme_name": result["records"][0].get("enzyme", {}).get("name"),
                        "data_points": result["count"],
                        "kcat_avg": sum(kcat_values) / len(kcat_values) if kcat_values else None,
                        "kcat_range": [min(kcat_values), max(kcat_values)] if kcat_values else None,
                        "km_avg": sum(km_values) / len(km_values) if km_values else None,
                        "km_range": [min(km_values), max(km_values)] if km_values else None,
                    })

            return json.dumps({"comparison": comparison}, indent=2)

        except Exception as e:
            logger.error(f"Error comparing enzymes: {e}")
            return json.dumps({"error": str(e)})

    def get_statistics(self) -> str:
        """
        Get overall knowledge base statistics.

        Returns:
            JSON string with statistics about the knowledge base

        Example:
            result = get_statistics()
        """
        self._ensure_loaded()

        try:
            stats = self.api.get_statistics()
            return json.dumps(stats, indent=2)

        except Exception as e:
            logger.error(f"Error getting statistics: {e}")
            return json.dumps({"error": str(e)})

    def get_tool_schemas(self) -> List[Dict[str, Any]]:
        """
        Get OpenAI-compatible function schemas for all tools.

        Returns:
            List of function schemas for Autogen registration
        """
        return [
            {
                "type": "function",
                "function": {
                    "name": "search_enzyme_by_ec",
                    "description": "Search for enzyme information by EC number. Returns enzyme names, organisms, reactions, and kinetic parameters.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "ec_number": {
                                "type": "string",
                                "description": "EC number in format X.X.X.X (e.g., '1.1.1.1')",
                            },
                            "include_sources": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "Optional list of data sources to filter (e.g., ['BRENDA', 'BKMS'])",
                            },
                        },
                        "required": ["ec_number"],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "search_reactions_by_compound",
                    "description": "Search for all enzymatic reactions involving a specific compound as substrate or product.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "compound_name": {
                                "type": "string",
                                "description": "Name of the compound (e.g., 'glucose', 'ATP')",
                            },
                            "role": {
                                "type": "string",
                                "enum": ["substrate", "product"],
                                "description": "Filter by role: 'substrate' or 'product'",
                            },
                        },
                        "required": ["compound_name"],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "find_retrosynthesis_pathway",
                    "description": "Find enzymatic retrosynthesis pathways for a target compound. Identifies enzymes that can produce the target from available substrates.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "target_compound": {
                                "type": "string",
                                "description": "Target compound to synthesize",
                            },
                            "available_substrates": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of available starting materials",
                            },
                        },
                        "required": ["target_compound"],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "get_kinetic_parameters",
                    "description": "Retrieve kinetic parameters (kcat, Km) for enzymes with optional filters for EC number, substrate, or organism.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "ec_number": {
                                "type": "string",
                                "description": "EC number filter",
                            },
                            "substrate": {
                                "type": "string",
                                "description": "Substrate name filter",
                            },
                            "organism": {
                                "type": "string",
                                "description": "Organism name filter",
                            },
                        },
                        "required": [],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "get_enzyme_sequence",
                    "description": "Retrieve amino acid sequences for enzymes by EC number, useful for protein engineering.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "ec_number": {
                                "type": "string",
                                "description": "EC number",
                            },
                            "organism": {
                                "type": "string",
                                "description": "Optional organism filter",
                            },
                        },
                        "required": ["ec_number"],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "compare_enzymes",
                    "description": "Compare kinetic parameters of multiple enzymes to select the best candidate.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "ec_numbers": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of EC numbers to compare",
                            },
                        },
                        "required": ["ec_numbers"],
                    },
                },
            },
            {
                "type": "function",
                "function": {
                    "name": "get_statistics",
                    "description": "Get overall statistics about the knowledge base including total records, EC numbers, organisms, etc.",
                    "parameters": {
                        "type": "object",
                        "properties": {},
                        "required": [],
                    },
                },
            },
        ]
