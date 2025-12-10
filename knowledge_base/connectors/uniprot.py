"""UniProt protein database connector (REST API)."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

import requests

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class UniProtConnector(BaseConnector):
    """Connector for UniProt REST API."""

    SOURCE_NAME = "UniProt"
    ENTITY_TYPE = "protein"
    BASE_URL = "https://rest.uniprot.org"

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        rate_limit: float = 1.0,
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, rate_limit=rate_limit, **kwargs)
        self._session: Optional[requests.Session] = None

    def connect(self) -> bool:
        """Initialize HTTP session."""
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": "ChemoenzymaticKnowledgeBase/1.0",
            "Accept": "application/json",
        })
        return True

    def disconnect(self) -> None:
        """Close HTTP session."""
        if self._session:
            self._session.close()
            self._session = None

    def is_available(self) -> bool:
        """Check if UniProt API is accessible."""
        if not self._session:
            return False
        try:
            resp = self._session.get(f"{self.BASE_URL}/uniprotkb/P00533", timeout=10)
            return resp.status_code in (200, 301, 302)
        except Exception:
            return False

    def _get(self, endpoint: str, params: Optional[Dict] = None) -> Optional[Dict]:
        """Make GET request to UniProt API."""
        if not self._session:
            return None

        self._rate_limit_wait()

        try:
            url = f"{self.BASE_URL}/{endpoint}"
            resp = self._session.get(url, params=params, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            logger.warning(f"UniProt request failed: {resp.status_code}")
            return None
        except Exception as e:
            logger.error(f"UniProt request error: {e}")
            return None

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query proteins by EC number."""
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

        # Search by EC number
        data = self._get("uniprotkb/search", params={
            "query": f"ec:{ec_number}",
            "format": "json",
            "size": 100,
        })

        if not data or "results" not in data:
            return ConnectorResult(
                success=False,
                error=f"No UniProt entries found for EC {ec_number}",
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        records = [self.to_unified_schema(entry) for entry in data["results"]]
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
        """Query proteins by name or accession."""
        start = time.time()

        if id_type == "accession":
            # Direct lookup by accession
            data = self._get(f"uniprotkb/{identifier}")
            if data:
                return ConnectorResult(
                    success=True,
                    data=[self.to_unified_schema(data)],
                    source=self.SOURCE_NAME,
                    query_time=time.time() - start,
                )
            return ConnectorResult(
                success=False,
                error=f"Protein not found: {identifier}",
                source=self.SOURCE_NAME,
            )
        else:
            # Search by name/keyword
            data = self._get("uniprotkb/search", params={
                "query": identifier,
                "format": "json",
                "size": 50,
            })

            if data and "results" in data:
                records = [self.to_unified_schema(entry) for entry in data["results"]]
                return ConnectorResult(
                    success=True,
                    data=records,
                    source=self.SOURCE_NAME,
                    query_time=time.time() - start,
                )

            return ConnectorResult(
                success=False,
                error=f"No results for: {identifier}",
                source=self.SOURCE_NAME,
            )

    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> ConnectorResult:
        """Query by substrate (search annotations)."""
        start = time.time()
        all_terms = (substrates or []) + (products or [])

        if not all_terms:
            return ConnectorResult(success=True, data=[], source=self.SOURCE_NAME)

        # Search using compound names in annotations
        query = " AND ".join([f'annotation:("{term}")' for term in all_terms[:3]])
        data = self._get("uniprotkb/search", params={
            "query": query,
            "format": "json",
            "size": 50,
        })

        if data and "results" in data:
            records = [self.to_unified_schema(entry) for entry in data["results"]]
            return ConnectorResult(
                success=True,
                data=records,
                source=self.SOURCE_NAME,
                query_time=time.time() - start,
            )

        return ConnectorResult(success=True, data=[], source=self.SOURCE_NAME)

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over enzymes (not practical for full UniProt)."""
        # This would be too large - return empty iterator
        logger.warning("get_all_records not supported for UniProt (too large)")
        return iter([])

    def get_protein(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get protein by UniProt accession."""
        data = self._get(f"uniprotkb/{accession}")
        if data:
            return self.to_unified_schema(data)
        return None

    def get_sequence(self, accession: str) -> Optional[str]:
        """Get protein sequence by accession."""
        if not self._session:
            return None

        self._rate_limit_wait()

        try:
            url = f"{self.BASE_URL}/uniprotkb/{accession}.fasta"
            resp = self._session.get(url, timeout=30)
            if resp.status_code == 200:
                lines = resp.text.strip().split("\n")
                return "".join(lines[1:])  # Skip header
            return None
        except Exception as e:
            logger.error(f"Failed to get sequence: {e}")
            return None

    def batch_get_proteins(self, accessions: List[str]) -> List[Dict[str, Any]]:
        """Get multiple proteins by accession."""
        records = []
        for acc in accessions:
            protein = self.get_protein(acc)
            if protein:
                records.append(protein)
        return records

    def to_unified_schema(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """Convert UniProt entry to unified schema."""
        accession = record.get("primaryAccession", "")
        organism = record.get("organism", {})
        organism_name = organism.get("scientificName", "")

        # Extract EC numbers
        ec_numbers = []
        for comment in record.get("comments", []):
            if comment.get("commentType") == "CATALYTIC ACTIVITY":
                for ec in comment.get("reaction", {}).get("ecNumber", "").split():
                    if ec:
                        ec_numbers.append(ec)

        # Extract protein name
        protein_name = ""
        protein_desc = record.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        if rec_name:
            protein_name = rec_name.get("fullName", {}).get("value", "")

        # Extract PDB cross-references
        pdb_ids = []
        for xref in record.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                pdb_ids.append(xref.get("id", ""))

        # Get sequence
        sequence = record.get("sequence", {}).get("value", "")

        return {
            "id": f"UniProt:{accession}",
            "source": {
                "dataset": "UniProt",
                "raw_ids": {"uniprot_accession": accession},
            },
            "ec_numbers": ec_numbers,
            "primary_ec": ec_numbers[0] if ec_numbers else None,
            "enzyme": {
                "name": protein_name,
                "name_full": protein_name,
                "uniprot_ids": [accession],
                "sequence": sequence,
                "pdb_ids": pdb_ids,
                "organism": organism_name,
                "mutant": None,
                "mutation_flag": False,
            },
            "reaction": {
                "direction": "unknown",
                "equation_text": None,
                "substrates": [],
                "products": [],
                "cofactors": [],
            },
            "kinetics": {
                "kcat": {"value": None, "unit": None},
                "km": {"value": None, "unit": None},
                "kcat_over_km": {"value": None, "unit": None},
            },
            "conditions": {},
            "text": {
                "descriptor": protein_name,
            },
            "uniprot_specific": {
                "entry_type": record.get("entryType"),
                "gene_names": [g.get("geneName", {}).get("value") for g in record.get("genes", [])],
                "keywords": [kw.get("name") for kw in record.get("keywords", [])],
                "features": len(record.get("features", [])),
            },
        }

    def map_ids(
        self,
        from_db: str,
        to_db: str,
        ids: List[str],
    ) -> Dict[str, List[str]]:
        """Map IDs between databases using UniProt ID mapping."""
        if not self._session or not ids:
            return {}

        self._rate_limit_wait()

        try:
            # Submit mapping job
            resp = self._session.post(
                f"{self.BASE_URL}/idmapping/run",
                data={
                    "from": from_db,
                    "to": to_db,
                    "ids": ",".join(ids),
                }
            )

            if resp.status_code != 200:
                return {}

            job_id = resp.json().get("jobId")
            if not job_id:
                return {}

            # Poll for results
            for _ in range(30):
                time.sleep(1)
                status_resp = self._session.get(
                    f"{self.BASE_URL}/idmapping/status/{job_id}"
                )
                if status_resp.status_code == 200:
                    status = status_resp.json()
                    if "results" in status or status.get("jobStatus") == "FINISHED":
                        break

            # Get results
            results_resp = self._session.get(
                f"{self.BASE_URL}/idmapping/results/{job_id}"
            )

            if results_resp.status_code == 200:
                data = results_resp.json()
                mapping = {}
                for result in data.get("results", []):
                    from_id = result.get("from", "")
                    to_id = result.get("to", {})
                    if isinstance(to_id, dict):
                        to_id = to_id.get("primaryAccession", "")
                    if from_id not in mapping:
                        mapping[from_id] = []
                    mapping[from_id].append(to_id)
                return mapping

        except Exception as e:
            logger.error(f"ID mapping failed: {e}")

        return {}
