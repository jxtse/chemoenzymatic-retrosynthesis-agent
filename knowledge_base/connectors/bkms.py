"""BKMS (Biochemical Kinetic Model System) connector."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import pandas as pd

from .base import BaseConnector, ConnectorResult
from ..schema import UnifiedRecord

logger = logging.getLogger(__name__)


class BKMSConnector(BaseConnector):
    """Connector for BKMS reaction database (local TSV file)."""

    SOURCE_NAME = "BKMS"
    ENTITY_TYPE = "reaction"

    def __init__(
        self,
        file_path: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, **kwargs)
        self.file_path = file_path
        self._df: Optional[pd.DataFrame] = None
        self._ec_index: Dict[str, List[int]] = {}

    def connect(self) -> bool:
        """Load the BKMS data file."""
        if self.file_path is None or not self.file_path.exists():
            logger.error(f"BKMS file not found: {self.file_path}")
            return False

        try:
            self._df = pd.read_csv(self.file_path, sep="\t", low_memory=False)
            self._build_ec_index()
            logger.info(f"Loaded BKMS with {len(self._df)} reactions")
            return True
        except Exception as e:
            logger.error(f"Failed to load BKMS: {e}")
            return False

    def _build_ec_index(self) -> None:
        """Build index mapping EC numbers to row indices."""
        self._ec_index = {}
        if self._df is None:
            return

        for idx, ec_val in enumerate(self._df.get("EC_Number", [])):
            if pd.notna(ec_val):
                ecs = self._parse_ec_list(str(ec_val))
                for ec in ecs:
                    if ec not in self._ec_index:
                        self._ec_index[ec] = []
                    self._ec_index[ec].append(idx)

    @staticmethod
    def _parse_ec_list(ec_str: str) -> List[str]:
        """Parse EC number string (may contain multiple ECs)."""
        if not ec_str or ec_str.lower() == "nan":
            return []
        # Handle various formats: "1.1.1.1", "1.1.1.1, 2.2.2.2", etc.
        ecs = []
        for part in ec_str.replace(";", ",").split(","):
            ec = part.strip()
            if ec and ec.lower() != "nan":
                ecs.append(ec)
        return ecs

    def disconnect(self) -> None:
        """Clear loaded data."""
        self._df = None
        self._ec_index = {}

    def is_available(self) -> bool:
        """Check if BKMS data is loaded."""
        return self._df is not None

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query reactions by EC number."""
        start = time.time()

        if not self.is_available():
            return ConnectorResult(
                success=False,
                error="BKMS data not loaded",
                source=self.SOURCE_NAME,
            )

        # Check cache
        cache_key = self._cache_key(f"ec:{ec_number}")
        cached = self._get_from_cache(cache_key)
        if cached is not None:
            return ConnectorResult(
                success=True,
                data=cached,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
                cached=True,
            )

        # Query index
        indices = self._ec_index.get(ec_number, [])
        records = []
        for idx in indices:
            row = self._df.iloc[idx].to_dict()
            records.append(self.to_unified_schema(row, idx))

        self._save_to_cache(cache_key, records)

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
            metadata={"count": len(records)},
        )

    def query_by_compound(
        self,
        identifier: str,
        id_type: str = "name",
    ) -> ConnectorResult:
        """Query reactions by compound name in equation."""
        start = time.time()

        if not self.is_available():
            return ConnectorResult(
                success=False,
                error="BKMS data not loaded",
                source=self.SOURCE_NAME,
            )

        identifier_lower = identifier.lower()
        records = []

        for idx, row in self._df.iterrows():
            equation = str(row.get("Reaction", "")).lower()
            if identifier_lower in equation:
                records.append(self.to_unified_schema(row.to_dict(), int(idx)))

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
            metadata={"count": len(records)},
        )

    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> ConnectorResult:
        """Query reactions by substrates/products."""
        start = time.time()

        if not self.is_available():
            return ConnectorResult(
                success=False,
                error="BKMS data not loaded",
                source=self.SOURCE_NAME,
            )

        substrates = substrates or []
        products = products or []
        records = []

        for idx, row in self._df.iterrows():
            equation = str(row.get("Reaction", ""))

            # Simple substring matching for substrates and products
            match = True
            for sub in substrates:
                if sub.lower() not in equation.lower():
                    match = False
                    break
            for prod in products:
                if prod.lower() not in equation.lower():
                    match = False
                    break

            if match:
                records.append(self.to_unified_schema(row.to_dict(), int(idx)))

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
            metadata={"count": len(records)},
        )

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over all BKMS reactions."""
        if not self.is_available():
            return

        for idx, row in self._df.iterrows():
            yield self.to_unified_schema(row.to_dict(), int(idx))

    def to_unified_schema(
        self,
        record: Dict[str, Any],
        row_index: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Convert BKMS record to unified schema."""
        row_index = row_index or 0

        # Parse EC numbers
        ecs = self._parse_ec_list(str(record.get("EC_Number", "")))

        # Parse reaction equation
        equation = record.get("Reaction", "")
        reaction_parts = self._parse_equation(equation)

        return {
            "id": f"BKMS:{row_index}",
            "source": {
                "dataset": "BKMS",
                "row_index": row_index,
                "raw_ids": {
                    "bkms_reaction_id": record.get("Reaction_ID_BRENDA"),
                    "kegg_id": record.get("Reaction_ID_KEGG"),
                    "metacyc_id": record.get("Reaction_ID_MetaCyc"),
                    "sabiork_id": record.get("Reaction_ID_SABIO_RK"),
                },
            },
            "ec_numbers": ecs,
            "primary_ec": ecs[0] if ecs else None,
            "enzyme": {
                "name": record.get("Recommended_Name"),
                "name_full": None,
                "uniprot_ids": [],
                "sequence": None,
                "pdb_ids": [],
                "organism": None,
                "mutant": None,
                "mutation_flag": False,
            },
            "reaction": {
                "direction": reaction_parts["direction"],
                "equation_text": equation,
                "substrates": reaction_parts["substrates"],
                "products": reaction_parts["products"],
                "cofactors": [],
                "reaction_smarts": None,
                "mechanism_type": None,
            },
            "kinetics": {
                "kcat": {"value": None, "unit": None, "raw_text": None, "condition_id": None},
                "km": {"value": None, "unit": None, "raw_text": None, "condition_id": None},
                "kcat_over_km": {"value": None, "unit": None},
            },
            "conditions": {
                "id": f"cond_BKMS_{row_index}",
                "pH": None,
                "temperature": None,
                "temperature_unit": None,
                "ionic_strength": None,
                "buffer": None,
                "comments": None,
            },
            "text": {
                "descriptor": record.get("Remark"),
                "canonical_description": None,
                "notes": record.get("Commentary_KEGG") or record.get("Commentary_MetaCyc"),
            },
            "pathway": {
                "kegg_pathway": record.get("KEGG_Pathway_Name"),
                "metacyc_pathway": record.get("MetaCyc_Pathway_Name"),
            },
        }

    @staticmethod
    def _parse_equation(equation: str) -> Dict[str, Any]:
        """Parse BKMS reaction equation."""
        if not equation:
            return {"direction": "unknown", "substrates": [], "products": []}

        direction = "unknown"
        separators = ["<=>", "=>", "<=", "->", "="]
        left, right = equation, ""

        for sep in separators:
            if sep in equation:
                parts = equation.split(sep, 1)
                left, right = parts[0], parts[1] if len(parts) > 1 else ""
                if sep in ("=>", "->", "="):
                    direction = "forward"
                elif sep == "<=":
                    direction = "reverse"
                else:
                    direction = "reversible"
                break

        def parse_side(side: str, role: str) -> List[Dict[str, Any]]:
            compounds = []
            for part in side.split("+"):
                part = part.strip()
                if not part:
                    continue
                stoich = 1.0
                tokens = part.split()
                if tokens and tokens[0].replace(".", "", 1).isdigit():
                    try:
                        stoich = float(tokens[0])
                        part = " ".join(tokens[1:]) or tokens[0]
                    except ValueError:
                        pass
                compounds.append({
                    "name": part,
                    "name_full": None,
                    "role": role,
                    "stoichiometry": stoich,
                    "smiles": None,
                    "cid": None,
                    "brenda_id": None,
                })
            return compounds

        return {
            "direction": direction,
            "substrates": parse_side(left, "substrate"),
            "products": parse_side(right, "product"),
        }
