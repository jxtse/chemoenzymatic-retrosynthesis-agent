"""KEGG database connector (REST API)."""

from __future__ import annotations

import logging
import re
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import requests

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class KEGGConnector(BaseConnector):
    """Connector for KEGG REST API."""

    SOURCE_NAME = "KEGG"
    ENTITY_TYPE = "pathway"
    BASE_URL = "https://rest.kegg.jp"

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        rate_limit: float = 0.2,  # KEGG requests 5 sec between requests
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, rate_limit=rate_limit, **kwargs)
        self._session: Optional[requests.Session] = None

    def connect(self) -> bool:
        """Initialize HTTP session."""
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": "ChemoenzymaticKnowledgeBase/1.0"
        })
        return True

    def disconnect(self) -> None:
        """Close HTTP session."""
        if self._session:
            self._session.close()
            self._session = None

    def is_available(self) -> bool:
        """Check if KEGG API is accessible."""
        if not self._session:
            return False
        try:
            resp = self._session.get(f"{self.BASE_URL}/info/kegg", timeout=10)
            return resp.status_code == 200
        except Exception:
            return False

    def _get(self, endpoint: str) -> Optional[str]:
        """Make GET request to KEGG API."""
        if not self._session:
            return None

        self._rate_limit_wait()

        try:
            url = f"{self.BASE_URL}/{endpoint}"
            resp = self._session.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.text
            logger.warning(f"KEGG request failed: {resp.status_code}")
            return None
        except Exception as e:
            logger.error(f"KEGG request error: {e}")
            return None

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query KEGG enzyme by EC number."""
        start = time.time()

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

        # Query enzyme entry
        text = self._get(f"get/ec:{ec_number}")
        if not text:
            return ConnectorResult(
                success=False,
                error=f"No KEGG entry found for EC {ec_number}",
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        record = self._parse_enzyme_entry(text, ec_number)

        # Get linked reactions
        reactions_text = self._get(f"link/reaction/ec:{ec_number}")
        if reactions_text:
            record["linked_reactions"] = self._parse_link_result(reactions_text)

        # Get linked pathways
        pathways_text = self._get(f"link/pathway/ec:{ec_number}")
        if pathways_text:
            record["linked_pathways"] = self._parse_link_result(pathways_text)

        self._save_to_cache(cache_key, record)

        return ConnectorResult(
            success=True,
            data=record,
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
        )

    def query_by_compound(
        self,
        identifier: str,
        id_type: str = "name",
    ) -> ConnectorResult:
        """Query KEGG compound."""
        start = time.time()

        if id_type == "name":
            # Search by name
            text = self._get(f"find/compound/{identifier}")
            if not text:
                return ConnectorResult(
                    success=False,
                    error=f"No compound found: {identifier}",
                    source=self.SOURCE_NAME,
                    query_time=time.time() - start,
                )

            compounds = []
            for line in text.strip().split("\n"):
                if line:
                    parts = line.split("\t", 1)
                    cpd_id = parts[0] if parts else ""
                    names = parts[1] if len(parts) > 1 else ""
                    compounds.append({
                        "kegg_id": cpd_id,
                        "names": names.split("; "),
                    })
            return ConnectorResult(
                success=True,
                data=compounds,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )
        elif id_type == "kegg":
            # Get by KEGG ID
            text = self._get(f"get/{identifier}")
            if not text:
                return ConnectorResult(
                    success=False,
                    error=f"No entry found: {identifier}",
                    source=self.SOURCE_NAME,
                    query_time=time.time() - start,
                )
            compound = self._parse_compound_entry(text, identifier)
            return ConnectorResult(
                success=True,
                data=compound,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        return ConnectorResult(
            success=False,
            error=f"Unsupported id_type: {id_type}",
            source=self.SOURCE_NAME,
        )

    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> ConnectorResult:
        """Search for reactions by compounds."""
        start = time.time()

        # Find reactions containing compounds
        compound_ids = (substrates or []) + (products or [])
        if not compound_ids:
            return ConnectorResult(
                success=True,
                data=[],
                source=self.SOURCE_NAME,
            )

        reaction_sets = []
        for cpd in compound_ids:
            text = self._get(f"link/reaction/{cpd}")
            if text:
                reactions = set(self._parse_link_result(text))
                reaction_sets.append(reactions)

        # Find intersection
        if reaction_sets:
            common_reactions = reaction_sets[0]
            for rs in reaction_sets[1:]:
                common_reactions &= rs

            # Get reaction details
            reactions = []
            for rxn_id in list(common_reactions)[:20]:  # Limit to 20
                rxn_text = self._get(f"get/{rxn_id}")
                if rxn_text:
                    reactions.append(self._parse_reaction_entry(rxn_text, rxn_id))

            return ConnectorResult(
                success=True,
                data=reactions,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        return ConnectorResult(
            success=True,
            data=[],
            source=self.SOURCE_NAME,
            query_time=time.time() - start,
        )

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over all KEGG enzymes (EC list)."""
        # Get all EC numbers
        text = self._get("list/enzyme")
        if not text:
            return

        for line in text.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t", 1)
            ec_id = parts[0].replace("ec:", "")

            result = self.query_by_ec(ec_id)
            if result.success and result.data:
                yield result.data

    def _parse_enzyme_entry(self, text: str, ec_number: str) -> Dict[str, Any]:
        """Parse KEGG enzyme flat file format."""
        entry = {
            "id": f"KEGG:{ec_number}",
            "source": {"dataset": "KEGG", "ec_number": ec_number},
            "ec_numbers": [ec_number],
            "primary_ec": ec_number,
            "name": None,
            "names": [],
            "reaction_text": None,
            "substrates": [],
            "products": [],
            "genes": [],
            "orthology": [],
            "dblinks": {},
        }

        current_field = None
        current_value = []

        for line in text.split("\n"):
            if not line:
                continue

            if line.startswith(" "):
                # Continuation of previous field
                if current_field:
                    current_value.append(line.strip())
            else:
                # Save previous field
                if current_field and current_value:
                    self._set_enzyme_field(entry, current_field, current_value)

                # New field
                parts = line.split(None, 1)
                current_field = parts[0] if parts else None
                current_value = [parts[1].strip()] if len(parts) > 1 else []

        # Save last field
        if current_field and current_value:
            self._set_enzyme_field(entry, current_field, current_value)

        return entry

    def _set_enzyme_field(self, entry: Dict, field: str, values: List[str]) -> None:
        """Set enzyme entry field."""
        combined = " ".join(values)

        if field == "NAME":
            entry["names"] = [n.strip().rstrip(";") for n in combined.split(";") if n.strip()]
            entry["name"] = entry["names"][0] if entry["names"] else None
        elif field == "SYSNAME":
            entry["systematic_name"] = combined
        elif field == "REACTION":
            entry["reaction_text"] = combined
        elif field == "SUBSTRATE":
            entry["substrates"] = [s.strip().rstrip(";") for s in combined.split(";") if s.strip()]
        elif field == "PRODUCT":
            entry["products"] = [p.strip().rstrip(";") for p in combined.split(";") if p.strip()]
        elif field == "GENES":
            entry["genes"].extend(values)
        elif field == "ORTHOLOGY":
            entry["orthology"].extend(values)
        elif field == "DBLINKS":
            for val in values:
                if ":" in val:
                    db, ids = val.split(":", 1)
                    entry["dblinks"][db.strip()] = ids.strip()

    def _parse_compound_entry(self, text: str, cpd_id: str) -> Dict[str, Any]:
        """Parse KEGG compound entry."""
        entry = {
            "id": f"KEGG:{cpd_id}",
            "kegg_id": cpd_id,
            "names": [],
            "formula": None,
            "mass": None,
            "dblinks": {},
        }

        current_field = None
        current_value = []

        for line in text.split("\n"):
            if not line:
                continue

            if line.startswith(" "):
                if current_field:
                    current_value.append(line.strip())
            else:
                if current_field and current_value:
                    combined = " ".join(current_value)
                    if current_field == "NAME":
                        entry["names"] = [n.strip().rstrip(";") for n in combined.split(";")]
                    elif current_field == "FORMULA":
                        entry["formula"] = combined
                    elif current_field == "EXACT_MASS":
                        try:
                            entry["mass"] = float(combined)
                        except ValueError:
                            pass
                    elif current_field == "DBLINKS":
                        for val in current_value:
                            if ":" in val:
                                db, ids = val.split(":", 1)
                                entry["dblinks"][db.strip()] = ids.strip()

                parts = line.split(None, 1)
                current_field = parts[0] if parts else None
                current_value = [parts[1].strip()] if len(parts) > 1 else []

        return entry

    def _parse_reaction_entry(self, text: str, rxn_id: str) -> Dict[str, Any]:
        """Parse KEGG reaction entry."""
        entry = {
            "id": f"KEGG:{rxn_id}",
            "kegg_id": rxn_id,
            "name": None,
            "equation": None,
            "ec_numbers": [],
            "enzymes": [],
        }

        for line in text.split("\n"):
            if line.startswith("NAME"):
                entry["name"] = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else None
            elif line.startswith("DEFINITION"):
                entry["equation"] = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else None
            elif line.startswith("EQUATION"):
                entry["equation_raw"] = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else None
            elif line.startswith("ENZYME"):
                ecs = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
                entry["ec_numbers"] = [e.strip() for e in ecs.split() if e.strip()]

        return entry

    def _parse_link_result(self, text: str) -> List[str]:
        """Parse KEGG link query result."""
        results = []
        for line in text.strip().split("\n"):
            if "\t" in line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    results.append(parts[1])
        return results

    def to_unified_schema(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """Convert KEGG record to unified schema."""
        ec = record.get("primary_ec", "")

        return {
            "id": record.get("id", f"KEGG:{ec}"),
            "source": {
                "dataset": "KEGG",
                "raw_ids": {"kegg_ec": ec},
            },
            "ec_numbers": record.get("ec_numbers", []),
            "primary_ec": ec,
            "enzyme": {
                "name": record.get("name"),
                "name_full": record.get("systematic_name"),
                "uniprot_ids": [],
                "sequence": None,
                "pdb_ids": [],
                "organism": None,
                "mutant": None,
                "mutation_flag": False,
            },
            "reaction": {
                "direction": "unknown",
                "equation_text": record.get("reaction_text"),
                "substrates": [{"name": s, "role": "substrate"} for s in record.get("substrates", [])],
                "products": [{"name": p, "role": "product"} for p in record.get("products", [])],
                "cofactors": [],
            },
            "kinetics": {
                "kcat": {"value": None, "unit": None},
                "km": {"value": None, "unit": None},
                "kcat_over_km": {"value": None, "unit": None},
            },
            "conditions": {},
            "text": {"descriptor": None},
            "kegg_specific": {
                "linked_reactions": record.get("linked_reactions", []),
                "linked_pathways": record.get("linked_pathways", []),
                "dblinks": record.get("dblinks", {}),
                "genes": record.get("genes", []),
            },
        }

    # Utility methods for common KEGG operations
    def get_pathway(self, pathway_id: str) -> Optional[Dict[str, Any]]:
        """Get pathway details."""
        text = self._get(f"get/{pathway_id}")
        if not text:
            return None
        # Parse pathway entry
        return {"id": pathway_id, "raw_text": text}

    def convert_ids(self, source_db: str, target_db: str, ids: List[str]) -> Dict[str, str]:
        """Convert IDs between databases."""
        mapping = {}
        for id_val in ids:
            text = self._get(f"conv/{target_db}/{source_db}:{id_val}")
            if text:
                for line in text.strip().split("\n"):
                    if "\t" in line:
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            mapping[id_val] = parts[1]
        return mapping
