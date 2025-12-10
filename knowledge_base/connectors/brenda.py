"""BRENDA enzyme database connector (local parquet + SOAP API)."""

from __future__ import annotations

import hashlib
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import pandas as pd

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class BRENDAConnector(BaseConnector):
    """Connector for BRENDA enzyme database."""

    SOURCE_NAME = "BRENDA"
    ENTITY_TYPE = "kinetics"

    def __init__(
        self,
        file_path: Optional[Path] = None,
        api_url: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        cache_dir: Optional[Path] = None,
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, **kwargs)
        self.file_path = file_path
        self.api_url = api_url or "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
        self.email = email or os.environ.get("BRENDA_EMAIL")
        self.password = password or os.environ.get("BRENDA_PASSWORD")

        self._df: Optional[pd.DataFrame] = None
        self._ec_index: Dict[str, List[int]] = {}
        self._client = None

    def connect(self) -> bool:
        """Load local file and optionally connect to SOAP API."""
        success = False

        # Load local parquet file
        if self.file_path and self.file_path.exists():
            try:
                self._df = pd.read_parquet(self.file_path)
                self._build_ec_index()
                logger.info(f"Loaded BRENDA parquet with {len(self._df)} records")
                success = True
            except Exception as e:
                logger.error(f"Failed to load BRENDA parquet: {e}")

        # Connect to SOAP API if credentials available
        if self.email and self.password:
            try:
                from zeep import Client
                self._client = Client(self.api_url)
                logger.info("Connected to BRENDA SOAP API")
                success = True
            except Exception as e:
                logger.warning(f"Failed to connect to BRENDA API: {e}")

        return success

    def _build_ec_index(self) -> None:
        """Build index mapping EC numbers to row indices."""
        self._ec_index = {}
        if self._df is None:
            return

        for idx, ec_val in enumerate(self._df.get("ec", [])):
            if pd.notna(ec_val):
                ec = str(ec_val).strip()
                if ec not in self._ec_index:
                    self._ec_index[ec] = []
                self._ec_index[ec].append(idx)

    def disconnect(self) -> None:
        """Clear loaded data."""
        self._df = None
        self._ec_index = {}
        self._client = None

    def is_available(self) -> bool:
        """Check if BRENDA data is loaded or API is connected."""
        return self._df is not None or self._client is not None

    def _get_credentials_hash(self) -> str:
        """Get hash of credentials for API calls."""
        if not self.email or not self.password:
            return ""
        cred_str = f"{self.email},{self.password}"
        return hashlib.sha256(cred_str.encode()).hexdigest()

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query kinetics by EC number."""
        start = time.time()
        records = []

        # Query local data first
        if self._df is not None:
            indices = self._ec_index.get(ec_number, [])
            for idx in indices:
                row = self._df.iloc[idx].to_dict()
                records.append(self.to_unified_schema(row, idx))

        # Query API if available and local data is empty
        if not records and self._client and self.email and self.password:
            try:
                api_records = self._query_api_km(ec_number)
                records.extend(api_records)
            except Exception as e:
                logger.warning(f"BRENDA API query failed: {e}")

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
            metadata={"count": len(records)},
        )

    def _query_api_km(self, ec_number: str) -> List[Dict[str, Any]]:
        """Query BRENDA SOAP API for Km values."""
        if not self._client:
            return []

        self._rate_limit_wait()

        try:
            cred_hash = self._get_credentials_hash()
            result = self._client.service.getKmValue(
                self.email,
                self.password,
                f"ecNumber*{ec_number}",
            )

            records = []
            for item in result or []:
                records.append({
                    "id": f"BRENDA_API:{ec_number}:{len(records)}",
                    "source": {"dataset": "BRENDA_API", "row_index": len(records)},
                    "ec_numbers": [ec_number],
                    "primary_ec": ec_number,
                    "enzyme": {
                        "name": None,
                        "organism": getattr(item, "organism", None),
                    },
                    "kinetics": {
                        "km": {
                            "value": getattr(item, "kmValue", None),
                            "unit": "mM",
                        },
                    },
                    "reaction": {
                        "substrates": [{"name": getattr(item, "substrate", None)}] if getattr(item, "substrate", None) else [],
                    },
                })
            return records
        except Exception as e:
            logger.error(f"BRENDA API error: {e}")
            return []

    def query_by_compound(
        self,
        identifier: str,
        id_type: str = "name",
    ) -> ConnectorResult:
        """Query kinetics by substrate name."""
        start = time.time()
        records = []

        if self._df is not None:
            identifier_lower = identifier.lower()
            for idx, row in self._df.iterrows():
                substrate = str(row.get("substrate", "")).lower()
                if identifier_lower in substrate:
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
        """Query by substrate (BRENDA doesn't have full reactions)."""
        if substrates:
            # BRENDA mainly has substrate-specific data
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
        """Iterate over all BRENDA records."""
        if self._df is None:
            return

        for idx, row in self._df.iterrows():
            yield self.to_unified_schema(row.to_dict(), int(idx))

    def to_unified_schema(
        self,
        record: Dict[str, Any],
        row_index: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Convert BRENDA record to unified schema."""
        row_index = row_index or 0

        ec = str(record.get("ec", "")).strip()
        ecs = [ec] if ec and ec.lower() != "nan" else []

        # Handle UniProt accessions
        accessions = record.get("accessions", "")
        uniprot_ids = []
        if pd.notna(accessions) and accessions:
            if isinstance(accessions, str):
                uniprot_ids = [a.strip() for a in accessions.split(",") if a.strip()]
            elif isinstance(accessions, list):
                uniprot_ids = accessions

        # Parse kinetic values
        kcat_val = record.get("turnover_number")
        km_val = record.get("km_value")
        kcat_km_val = record.get("kcat_km")

        return {
            "id": f"BRENDA:{row_index}",
            "source": {
                "dataset": "BRENDA",
                "row_index": row_index,
                "raw_ids": {
                    "pmid": record.get("pmid"),
                    "protein_id": record.get("protein_id"),
                },
            },
            "ec_numbers": ecs,
            "primary_ec": ecs[0] if ecs else None,
            "enzyme": {
                "name": record.get("enzyme"),
                "name_full": None,
                "uniprot_ids": uniprot_ids,
                "sequence": None,
                "pdb_ids": [],
                "organism": record.get("organism_name"),
                "mutant": record.get("mutant"),
                "mutation_flag": bool(record.get("mutant")) and str(record.get("mutant")).upper() != "WT",
            },
            "reaction": {
                "direction": "unknown",
                "equation_text": None,
                "substrates": [{
                    "name": record.get("substrate"),
                    "name_full": None,
                    "role": "substrate",
                    "stoichiometry": 1.0,
                    "smiles": None,
                    "cid": None,
                    "brenda_id": None,
                }] if record.get("substrate") else [],
                "products": [],
                "cofactors": [],
                "reaction_smarts": None,
                "mechanism_type": None,
            },
            "kinetics": {
                "kcat": {
                    "value": float(kcat_val) if pd.notna(kcat_val) else None,
                    "unit": "s^-1",
                    "raw_text": str(kcat_val) if pd.notna(kcat_val) else None,
                    "condition_id": f"cond_BRENDA_{row_index}",
                },
                "km": {
                    "value": float(km_val) if pd.notna(km_val) else None,
                    "unit": "mM",
                    "raw_text": str(km_val) if pd.notna(km_val) else None,
                    "condition_id": f"cond_BRENDA_{row_index}",
                },
                "kcat_over_km": {
                    "value": float(kcat_km_val) if pd.notna(kcat_km_val) else None,
                    "unit": "s^-1 mM^-1",
                },
            },
            "conditions": {
                "id": f"cond_BRENDA_{row_index}",
                "pH": float(record.get("pH")) if pd.notna(record.get("pH")) else None,
                "temperature": float(record.get("temperature")) if pd.notna(record.get("temperature")) else None,
                "temperature_unit": "C",
                "ionic_strength": None,
                "buffer": None,
                "comments": record.get("comments"),
            },
            "text": {
                "descriptor": record.get("comments"),
                "canonical_description": None,
                "notes": record.get("organism_comments"),
            },
        }

    def get_kcat_values(
        self,
        ec_number: str,
        organism: str = "*",
        substrate: str = "*",
    ) -> List[Dict[str, Any]]:
        """Get kcat values for an enzyme (API method)."""
        if not self._client or not self.email:
            return []

        self._rate_limit_wait()

        try:
            query = f"ecNumber*{ec_number}"
            if organism != "*":
                query += f"#organism*{organism}"
            if substrate != "*":
                query += f"#substrate*{substrate}"

            result = self._client.service.getTurnoverNumber(
                self.email, self.password, query
            )
            return [
                {
                    "ec": ec_number,
                    "organism": getattr(item, "organism", None),
                    "substrate": getattr(item, "substrate", None),
                    "kcat": getattr(item, "turnoverNumber", None),
                    "comment": getattr(item, "commentary", None),
                }
                for item in result or []
            ]
        except Exception as e:
            logger.error(f"Failed to get kcat values: {e}")
            return []
