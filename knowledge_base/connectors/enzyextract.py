"""EnzyExtract database connector (local parquet file)."""

from __future__ import annotations

import ast
import logging
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import pandas as pd

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class EnzyExtractConnector(BaseConnector):
    """Connector for EnzyExtract kinetics database."""

    SOURCE_NAME = "EnzyExtract"
    ENTITY_TYPE = "kinetics"

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
        """Load the EnzyExtract parquet file."""
        if self.file_path is None or not self.file_path.exists():
            logger.error(f"EnzyExtract file not found: {self.file_path}")
            return False

        try:
            self._df = pd.read_parquet(self.file_path)
            self._build_ec_index()
            logger.info(f"Loaded EnzyExtract with {len(self._df)} records")
            return True
        except Exception as e:
            logger.error(f"Failed to load EnzyExtract: {e}")
            return False

    def _build_ec_index(self) -> None:
        """Build index mapping EC numbers to row indices."""
        self._ec_index = {}
        if self._df is None:
            return

        for idx, ec_val in enumerate(self._df.get("enzyme_ecs", [])):
            ecs = self._parse_ec_list(ec_val)
            for ec in ecs:
                if ec not in self._ec_index:
                    self._ec_index[ec] = []
                self._ec_index[ec].append(idx)

    @staticmethod
    def _parse_ec_list(ec_val: Any) -> List[str]:
        """Parse EC number field (may be string repr of list)."""
        if pd.isna(ec_val):
            return []
        if isinstance(ec_val, list):
            return [str(e).strip() for e in ec_val if e and str(e).strip()]
        if isinstance(ec_val, str):
            text = ec_val.strip()
            if text.startswith("[") and text.endswith("]"):
                try:
                    parsed = ast.literal_eval(text)
                    if isinstance(parsed, list):
                        return [str(e).strip() for e in parsed if e]
                except (ValueError, SyntaxError):
                    pass
            # Fallback: split on commas/semicolons
            parts = [p.strip() for p in text.replace(";", ",").split(",")]
            return [p for p in parts if p and p.lower() != "nan"]
        return []

    def disconnect(self) -> None:
        """Clear loaded data."""
        self._df = None
        self._ec_index = {}

    def is_available(self) -> bool:
        """Check if data is loaded."""
        return self._df is not None

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query kinetics by EC number."""
        start = time.time()

        if not self.is_available():
            return ConnectorResult(
                success=False,
                error="EnzyExtract data not loaded",
                source=self.SOURCE_NAME,
            )

        indices = self._ec_index.get(ec_number, [])
        records = []
        for idx in indices:
            row = self._df.iloc[idx].to_dict()
            records.append(self.to_unified_schema(row, idx))

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
        """Query kinetics by substrate or compound identifier."""
        start = time.time()

        if not self.is_available():
            return ConnectorResult(
                success=False,
                error="EnzyExtract data not loaded",
                source=self.SOURCE_NAME,
            )

        records = []
        identifier_lower = identifier.lower()

        for idx, row in self._df.iterrows():
            match = False

            if id_type == "name":
                substrate = str(row.get("substrate", "")).lower()
                substrate_full = str(row.get("substrate_full", "")).lower()
                match = identifier_lower in substrate or identifier_lower in substrate_full
            elif id_type == "smiles":
                smiles = str(row.get("smiles", ""))
                match = identifier == smiles
            elif id_type == "cid":
                cid = str(row.get("cid", ""))
                match = identifier == cid

            if match:
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
        """Query by substrate (EnzyExtract has substrate-level data)."""
        if substrates:
            results = []
            for sub in substrates:
                result = self.query_by_compound(sub)
                if result.success and result.data:
                    results.extend(result.data)
            return ConnectorResult(
                success=True,
                data=results,
                source=self.SOURCE_NAME,
            )
        return ConnectorResult(
            success=True,
            data=[],
            source=self.SOURCE_NAME,
        )

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over all records."""
        if self._df is None:
            return

        for idx, row in self._df.iterrows():
            yield self.to_unified_schema(row.to_dict(), int(idx))

    def to_unified_schema(
        self,
        record: Dict[str, Any],
        row_index: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Convert EnzyExtract record to unified schema."""
        row_index = row_index or 0

        ecs = self._parse_ec_list(record.get("enzyme_ecs"))

        # Parse UniProt IDs
        uniprot_ids = self._parse_ec_list(record.get("uniprot"))
        pdb_ids = self._parse_ec_list(record.get("pdb"))

        # Parse cofactors
        cofactors = []
        cofactor_val = record.get("cofactors")
        if pd.notna(cofactor_val) and cofactor_val:
            cof_list = self._parse_ec_list(cofactor_val)
            for cof in cof_list:
                cofactors.append({
                    "name": cof,
                    "role": "cofactor",
                    "stoichiometry": 1.0,
                    "smiles": None,
                })

        # Parse kinetic values
        kcat_val = record.get("kcat_value")
        km_val = record.get("km_value")
        kcat_km_val = record.get("kcat_km")

        return {
            "id": f"EnzyExtract:{row_index}",
            "source": {
                "dataset": "EnzyExtract",
                "row_index": row_index,
                "raw_ids": {
                    "brenda_id": record.get("brenda_id"),
                    "cid": record.get("cid"),
                },
            },
            "ec_numbers": ecs,
            "primary_ec": ecs[0] if ecs else None,
            "enzyme": {
                "name": record.get("enzyme"),
                "name_full": record.get("enzyme_full"),
                "uniprot_ids": uniprot_ids,
                "sequence": record.get("sequence"),
                "pdb_ids": pdb_ids,
                "organism": record.get("organism"),
                "mutant": record.get("mutant"),
                "mutation_flag": bool(record.get("mutant")) and str(record.get("mutant")).upper() != "WT",
            },
            "reaction": {
                "direction": "unknown",
                "equation_text": None,
                "substrates": [{
                    "name": record.get("substrate"),
                    "name_full": record.get("substrate_full"),
                    "role": "substrate",
                    "stoichiometry": 1.0,
                    "smiles": record.get("smiles"),
                    "cid": record.get("cid"),
                    "brenda_id": record.get("brenda_id"),
                }] if record.get("substrate") else [],
                "products": [],
                "cofactors": cofactors,
                "reaction_smarts": None,
                "mechanism_type": None,
            },
            "kinetics": {
                "kcat": {
                    "value": float(kcat_val) if pd.notna(kcat_val) else None,
                    "unit": "s^-1",
                    "raw_text": record.get("kcat"),
                    "condition_id": f"cond_EnzyExtract_{row_index}",
                },
                "km": {
                    "value": float(km_val) if pd.notna(km_val) else None,
                    "unit": "mM",
                    "raw_text": record.get("km"),
                    "condition_id": f"cond_EnzyExtract_{row_index}",
                },
                "kcat_over_km": {
                    "value": float(kcat_km_val) if pd.notna(kcat_km_val) else None,
                    "unit": "s^-1 mM^-1",
                },
            },
            "conditions": {
                "id": f"cond_EnzyExtract_{row_index}",
                "pH": float(record.get("pH")) if pd.notna(record.get("pH")) else None,
                "temperature": float(record.get("temperature")) if pd.notna(record.get("temperature")) else None,
                "temperature_unit": "C",
                "ionic_strength": None,
                "buffer": None,
                "comments": record.get("solution"),
            },
            "text": {
                "descriptor": record.get("descriptor"),
                "canonical_description": record.get("canonical"),
                "notes": record.get("other"),
            },
        }
