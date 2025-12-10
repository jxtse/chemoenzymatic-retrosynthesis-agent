"""Unified JSON API for querying the knowledge base."""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class KnowledgeBaseAPI:
    """Unified API for querying the integrated knowledge base."""

    def __init__(self, kb_path: Path):
        """
        Initialize API with knowledge base file.

        Args:
            kb_path: Path to knowledge base JSONL or parquet file
        """
        self.kb_path = kb_path
        self._records: List[Dict[str, Any]] = []
        self._ec_index: Dict[str, List[int]] = {}
        self._compound_index: Dict[str, List[int]] = {}
        self._loaded = False

    def load(self) -> None:
        """Load knowledge base into memory and build indices."""
        if self._loaded:
            return

        logger.info(f"Loading knowledge base from {self.kb_path}...")

        if self.kb_path.suffix == ".jsonl" or self.kb_path.suffix == ".gz":
            self._load_jsonl()
        elif self.kb_path.suffix == ".parquet":
            self._load_parquet()
        else:
            raise ValueError(f"Unsupported file format: {self.kb_path.suffix}")

        self._build_indices()
        self._loaded = True

        logger.info(f"✓ Loaded {len(self._records)} records")
        logger.info(f"✓ Indexed {len(self._ec_index)} EC numbers")
        logger.info(f"✓ Indexed {len(self._compound_index)} compounds")

    def _load_jsonl(self) -> None:
        """Load JSONL knowledge base."""
        import gzip

        if self.kb_path.suffix == ".gz":
            opener = gzip.open
        else:
            opener = open

        with opener(self.kb_path, "rt", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    self._records.append(json.loads(line))

    def _load_parquet(self) -> None:
        """Load Parquet knowledge base."""
        import pandas as pd
        df = pd.read_parquet(self.kb_path)
        self._records = df.to_dict("records")

    def _build_indices(self) -> None:
        """Build search indices."""
        self._ec_index = defaultdict(list)
        self._compound_index = defaultdict(list)

        for idx, record in enumerate(self._records):
            # Index by EC number
            for ec in record.get("ec_numbers", []):
                self._ec_index[ec].append(idx)

            # Index by compound names
            for substrate in record.get("reaction", {}).get("substrates", []):
                name = substrate.get("name", "").lower()
                if name:
                    self._compound_index[name].append(idx)

            for product in record.get("reaction", {}).get("products", []):
                name = product.get("name", "").lower()
                if name:
                    self._compound_index[name].append(idx)

    def query_by_ec(
        self,
        ec_number: str,
        include_sources: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Query knowledge base by EC number.

        Args:
            ec_number: EC number (e.g., "1.1.1.1")
            include_sources: Filter by data sources (e.g., ["BKMS", "BRENDA"])

        Returns:
            Dict with query results
        """
        if not self._loaded:
            self.load()

        indices = self._ec_index.get(ec_number, [])
        records = [self._records[i] for i in indices]

        # Filter by sources if specified
        if include_sources:
            records = [
                r for r in records
                if r.get("source", {}).get("dataset") in include_sources
            ]

        return {
            "query": {"ec_number": ec_number},
            "count": len(records),
            "records": records,
        }

    def query_by_compound(
        self,
        compound_name: str,
        role: Optional[str] = None,  # "substrate", "product", or None for both
    ) -> Dict[str, Any]:
        """
        Query knowledge base by compound name.

        Args:
            compound_name: Compound name
            role: Filter by role (substrate/product)

        Returns:
            Dict with query results
        """
        if not self._loaded:
            self.load()

        compound_lower = compound_name.lower()
        indices = self._compound_index.get(compound_lower, [])

        # Also do fuzzy matching
        fuzzy_indices = set()
        for name, idx_list in self._compound_index.items():
            if compound_lower in name or name in compound_lower:
                fuzzy_indices.update(idx_list)

        all_indices = set(indices) | fuzzy_indices
        records = [self._records[i] for i in all_indices]

        # Filter by role if specified
        if role:
            filtered = []
            for record in records:
                reaction = record.get("reaction", {})
                items = reaction.get(f"{role}s", []) if role in ["substrate", "product"] else []
                for item in items:
                    if compound_lower in item.get("name", "").lower():
                        filtered.append(record)
                        break
            records = filtered

        return {
            "query": {"compound": compound_name, "role": role},
            "count": len(records),
            "records": records,
        }

    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Query by reaction (substrates and/or products).

        Args:
            substrates: List of substrate names
            products: List of product names

        Returns:
            Dict with query results
        """
        if not self._loaded:
            self.load()

        # Find records matching all substrates
        substrate_matches = None
        if substrates:
            for sub in substrates:
                result = self.query_by_compound(sub, role="substrate")
                record_ids = set(r["id"] for r in result["records"])
                if substrate_matches is None:
                    substrate_matches = record_ids
                else:
                    substrate_matches &= record_ids

        # Find records matching all products
        product_matches = None
        if products:
            for prod in products:
                result = self.query_by_compound(prod, role="product")
                record_ids = set(r["id"] for r in result["records"])
                if product_matches is None:
                    product_matches = record_ids
                else:
                    product_matches &= record_ids

        # Combine
        if substrate_matches is not None and product_matches is not None:
            matching_ids = substrate_matches & product_matches
        elif substrate_matches is not None:
            matching_ids = substrate_matches
        elif product_matches is not None:
            matching_ids = product_matches
        else:
            matching_ids = set()

        records = [r for r in self._records if r["id"] in matching_ids]

        return {
            "query": {"substrates": substrates, "products": products},
            "count": len(records),
            "records": records,
        }

    def get_kinetics(
        self,
        ec_number: Optional[str] = None,
        substrate: Optional[str] = None,
        organism: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Get kinetic parameters (kcat, Km) with optional filters.

        Args:
            ec_number: Filter by EC number
            substrate: Filter by substrate name
            organism: Filter by organism

        Returns:
            Dict with kinetic data
        """
        if not self._loaded:
            self.load()

        records = self._records

        # Filter by EC
        if ec_number:
            indices = self._ec_index.get(ec_number, [])
            records = [self._records[i] for i in indices]

        # Filter by substrate
        if substrate:
            substrate_lower = substrate.lower()
            records = [
                r for r in records
                if any(
                    substrate_lower in s.get("name", "").lower()
                    for s in r.get("reaction", {}).get("substrates", [])
                )
            ]

        # Filter by organism
        if organism:
            organism_lower = organism.lower()
            records = [
                r for r in records
                if organism_lower in r.get("enzyme", {}).get("organism", "").lower()
            ]

        # Extract kinetics
        kinetics = []
        for record in records:
            kin = record.get("kinetics", {})
            cond = record.get("conditions", {})

            kcat_val = kin.get("kcat", {}).get("value")
            km_val = kin.get("km", {}).get("value")

            if kcat_val is not None or km_val is not None:
                kinetics.append({
                    "id": record.get("id"),
                    "ec_number": record.get("primary_ec"),
                    "substrate": record.get("reaction", {}).get("substrates", [{}])[0].get("name"),
                    "organism": record.get("enzyme", {}).get("organism"),
                    "kcat": kcat_val,
                    "kcat_unit": kin.get("kcat", {}).get("unit"),
                    "km": km_val,
                    "km_unit": kin.get("km", {}).get("unit"),
                    "pH": cond.get("pH"),
                    "temperature": cond.get("temperature"),
                    "source": record.get("source", {}).get("dataset"),
                })

        return {
            "query": {
                "ec_number": ec_number,
                "substrate": substrate,
                "organism": organism,
            },
            "count": len(kinetics),
            "kinetics": kinetics,
        }

    def get_statistics(self) -> Dict[str, Any]:
        """Get knowledge base statistics."""
        if not self._loaded:
            self.load()

        # Count by source
        by_source = defaultdict(int)
        by_ec_class = defaultdict(int)
        organisms = set()
        has_kinetics = 0
        has_sequence = 0

        for record in self._records:
            source = record.get("source", {}).get("dataset", "Unknown")
            by_source[source] += 1

            # EC class (first digit)
            ec = record.get("primary_ec")
            if ec:
                ec_class = ec.split(".")[0] if "." in ec else ""
                if ec_class:
                    by_ec_class[f"EC {ec_class}"] += 1

            # Organism
            org = record.get("enzyme", {}).get("organism")
            if org:
                organisms.add(org)

            # Has kinetics
            kin = record.get("kinetics", {})
            if kin.get("kcat", {}).get("value") or kin.get("km", {}).get("value"):
                has_kinetics += 1

            # Has sequence
            if record.get("enzyme", {}).get("sequence"):
                has_sequence += 1

        return {
            "total_records": len(self._records),
            "unique_ec_numbers": len(self._ec_index),
            "unique_compounds": len(self._compound_index),
            "unique_organisms": len(organisms),
            "records_with_kinetics": has_kinetics,
            "records_with_sequence": has_sequence,
            "by_source": dict(by_source),
            "by_ec_class": dict(by_ec_class),
        }

    def export_to_json(
        self,
        output_path: Path,
        query_result: Dict[str, Any],
        indent: int = 2,
    ) -> None:
        """Export query result to JSON file."""
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(query_result, f, indent=indent, ensure_ascii=False)
        logger.info(f"✓ Exported to {output_path}")

    def __enter__(self) -> "KnowledgeBaseAPI":
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        pass
