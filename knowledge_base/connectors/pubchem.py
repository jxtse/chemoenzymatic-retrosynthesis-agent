"""PubChem compound database connector (REST API)."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import requests

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class PubChemConnector(BaseConnector):
    """Connector for PubChem PUG REST API."""

    SOURCE_NAME = "PubChem"
    ENTITY_TYPE = "compound"
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        rate_limit: float = 5.0,  # PubChem allows 5 req/sec
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, rate_limit=rate_limit, **kwargs)
        self._session: Optional[requests.Session] = None

    def connect(self) -> bool:
        """Initialize HTTP session."""
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": "ChemoenzymaticKnowledgeBase/1.0",
        })
        return True

    def disconnect(self) -> None:
        """Close HTTP session."""
        if self._session:
            self._session.close()
            self._session = None

    def is_available(self) -> bool:
        """Check if PubChem API is accessible."""
        if not self._session:
            return False
        try:
            resp = self._session.get(
                f"{self.BASE_URL}/compound/cid/2244/property/MolecularFormula/JSON",
                timeout=10
            )
            return resp.status_code == 200
        except Exception:
            return False

    def _get(self, endpoint: str) -> Optional[Dict]:
        """Make GET request to PubChem API."""
        if not self._session:
            return None

        self._rate_limit_wait()

        try:
            url = f"{self.BASE_URL}/{endpoint}"
            resp = self._session.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            logger.warning(f"PubChem request failed: {resp.status_code}")
            return None
        except Exception as e:
            logger.error(f"PubChem request error: {e}")
            return None

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query compounds associated with EC number (not directly supported)."""
        # PubChem doesn't directly support EC queries
        # Return empty result
        return ConnectorResult(
            success=True,
            data=[],
            source=self.SOURCE_NAME,
            metadata={"note": "EC queries not supported by PubChem"},
        )

    def query_by_compound(
        self,
        identifier: str,
        id_type: str = "name",
    ) -> ConnectorResult:
        """Query compound by various identifiers."""
        start = time.time()

        # Check cache
        cache_key = self._cache_key(f"{id_type}:{identifier}")
        cached = self._get_from_cache(cache_key)
        if cached is not None:
            return ConnectorResult(
                success=True,
                data=cached,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
                cached=True,
            )

        data = None
        if id_type == "name":
            data = self._get(f"compound/name/{identifier}/JSON")
        elif id_type == "smiles":
            data = self._get(f"compound/smiles/{identifier}/JSON")
        elif id_type == "inchi":
            data = self._get(f"compound/inchi/{identifier}/JSON")
        elif id_type == "cid":
            data = self._get(f"compound/cid/{identifier}/JSON")
        elif id_type == "inchikey":
            data = self._get(f"compound/inchikey/{identifier}/JSON")

        if not data:
            return ConnectorResult(
                success=False,
                error=f"No compound found: {identifier}",
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        # Get additional properties
        compounds = data.get("PC_Compounds", [])
        records = []

        for compound in compounds:
            cid = compound.get("id", {}).get("id", {}).get("cid")
            if cid:
                props = self._get_compound_properties(cid)
                record = self.to_unified_schema(compound, props)
                records.append(record)

        self._save_to_cache(cache_key, records)

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
            metadata={"count": len(records)},
        )

    def _get_compound_properties(self, cid: int) -> Dict[str, Any]:
        """Get compound properties by CID."""
        properties = [
            "MolecularFormula", "MolecularWeight", "CanonicalSMILES",
            "IsomericSMILES", "InChI", "InChIKey", "IUPACName",
            "XLogP", "TPSA", "HBondDonorCount", "HBondAcceptorCount",
            "RotatableBondCount", "HeavyAtomCount"
        ]
        prop_str = ",".join(properties)
        data = self._get(f"compound/cid/{cid}/property/{prop_str}/JSON")

        if data and "PropertyTable" in data:
            props = data["PropertyTable"].get("Properties", [])
            if props:
                return props[0]
        return {}

    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> ConnectorResult:
        """Query compounds by name (gather info for reaction components)."""
        start = time.time()
        all_compounds = (substrates or []) + (products or [])
        records = []

        for compound in all_compounds:
            result = self.query_by_compound(compound, "name")
            if result.success and result.data:
                records.extend(result.data)

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
        )

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Not practical for PubChem (100M+ compounds)."""
        logger.warning("get_all_records not supported for PubChem (too large)")
        return iter([])

    def similarity_search(
        self,
        smiles: str,
        threshold: float = 0.9,
        max_records: int = 100,
    ) -> List[Dict[str, Any]]:
        """Find similar compounds by SMILES."""
        data = self._get(
            f"compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON"
            f"?Threshold={int(threshold * 100)}&MaxRecords={max_records}"
        )

        if not data:
            return []

        cids = data.get("IdentifierList", {}).get("CID", [])
        records = []

        for cid in cids[:max_records]:
            result = self.query_by_compound(str(cid), "cid")
            if result.success and result.data:
                records.extend(result.data)

        return records

    def substructure_search(
        self,
        smiles: str,
        max_records: int = 100,
    ) -> List[Dict[str, Any]]:
        """Find compounds containing substructure."""
        data = self._get(
            f"compound/fastsubstructure/smiles/{smiles}/cids/JSON"
            f"?MaxRecords={max_records}"
        )

        if not data:
            return []

        cids = data.get("IdentifierList", {}).get("CID", [])
        records = []

        for cid in cids[:max_records]:
            result = self.query_by_compound(str(cid), "cid")
            if result.success and result.data:
                records.extend(result.data)

        return records

    def get_compound_bioactivity(self, cid: int) -> List[Dict[str, Any]]:
        """Get bioactivity data for a compound."""
        data = self._get(f"compound/cid/{cid}/assaysummary/JSON")

        if not data:
            return []

        table = data.get("Table", {})
        columns = table.get("Columns", {}).get("Column", [])
        rows = table.get("Row", [])

        bioactivities = []
        for row in rows:
            cells = row.get("Cell", [])
            if len(cells) >= len(columns):
                activity = dict(zip(columns, cells))
                bioactivities.append(activity)

        return bioactivities

    def to_unified_schema(
        self,
        compound: Dict[str, Any],
        properties: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Convert PubChem compound to unified schema."""
        properties = properties or {}

        cid = compound.get("id", {}).get("id", {}).get("cid", 0)

        # Extract SMILES from compound data
        smiles = properties.get("CanonicalSMILES", "")
        if not smiles:
            for prop in compound.get("props", []):
                if prop.get("urn", {}).get("label") == "SMILES":
                    if prop.get("urn", {}).get("name") == "Canonical":
                        smiles = prop.get("value", {}).get("sval", "")
                        break

        # Extract InChI
        inchi = properties.get("InChI", "")
        inchikey = properties.get("InChIKey", "")

        return {
            "id": f"PubChem:{cid}",
            "source": {
                "dataset": "PubChem",
                "raw_ids": {"cid": cid},
            },
            "compound": {
                "cid": cid,
                "name": properties.get("IUPACName", ""),
                "formula": properties.get("MolecularFormula", ""),
                "molecular_weight": properties.get("MolecularWeight"),
                "smiles": smiles,
                "isomeric_smiles": properties.get("IsomericSMILES", ""),
                "inchi": inchi,
                "inchikey": inchikey,
            },
            "properties": {
                "xlogp": properties.get("XLogP"),
                "tpsa": properties.get("TPSA"),
                "hbd": properties.get("HBondDonorCount"),
                "hba": properties.get("HBondAcceptorCount"),
                "rotatable_bonds": properties.get("RotatableBondCount"),
                "heavy_atoms": properties.get("HeavyAtomCount"),
            },
            # For compound records, enzyme/reaction fields are minimal
            "ec_numbers": [],
            "primary_ec": None,
            "enzyme": {},
            "reaction": {
                "substrates": [{
                    "name": properties.get("IUPACName", ""),
                    "smiles": smiles,
                    "cid": str(cid),
                }] if cid else [],
            },
            "kinetics": {},
            "conditions": {},
        }

    def get_synonyms(self, cid: int) -> List[str]:
        """Get compound synonyms."""
        data = self._get(f"compound/cid/{cid}/synonyms/JSON")
        if data:
            info = data.get("InformationList", {}).get("Information", [])
            if info:
                return info[0].get("Synonym", [])
        return []

    def smiles_to_cid(self, smiles: str) -> Optional[int]:
        """Convert SMILES to PubChem CID."""
        data = self._get(f"compound/smiles/{smiles}/cids/JSON")
        if data:
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0]
        return None

    def name_to_cid(self, name: str) -> Optional[int]:
        """Convert compound name to PubChem CID."""
        data = self._get(f"compound/name/{name}/cids/JSON")
        if data:
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0]
        return None
