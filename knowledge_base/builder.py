"""Knowledge base builder - integrates all data sources."""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from tqdm import tqdm

from .config import KnowledgeBaseConfig
from .connectors import (
    BKMSConnector,
    BRENDAConnector,
    EnzyExtractConnector,
    KEGGConnector,
    PubChemConnector,
    RetroBioCatConnector,
    UniProtConnector,
    USPTOConnector,
    ZINCConnector,
)
from .schema import SchemaValidator

logger = logging.getLogger(__name__)


class KnowledgeBaseBuilder:
    """Build unified knowledge base from multiple data sources."""

    def __init__(self, config: KnowledgeBaseConfig):
        self.config = config
        self.connectors: Dict[str, Any] = {}
        self.validator = SchemaValidator()

        # Statistics
        self.stats = {
            "total_records": 0,
            "by_source": defaultdict(int),
            "by_ec": defaultdict(int),
            "errors": [],
        }

    def initialize_connectors(self) -> None:
        """Initialize all enabled data source connectors."""
        logger.info("Initializing data source connectors...")

        # Local data sources
        if self.config.bkms.enabled and self.config.bkms.path:
            self.connectors["bkms"] = BKMSConnector(
                file_path=self.config.bkms.path,
                cache_dir=self.config.cache_dir / "bkms",
            )

        if self.config.brenda.enabled and self.config.brenda.path:
            self.connectors["brenda"] = BRENDAConnector(
                file_path=self.config.brenda.path,
                api_url=self.config.brenda.api_url,
                cache_dir=self.config.cache_dir / "brenda",
            )

        if self.config.enzyextract.enabled and self.config.enzyextract.path:
            self.connectors["enzyextract"] = EnzyExtractConnector(
                file_path=self.config.enzyextract.path,
                cache_dir=self.config.cache_dir / "enzyextract",
            )

        # API-based sources
        if self.config.kegg.enabled:
            self.connectors["kegg"] = KEGGConnector(
                cache_dir=self.config.cache_dir / "kegg",
                rate_limit=self.config.kegg.rate_limit,
            )

        if self.config.uniprot.enabled:
            self.connectors["uniprot"] = UniProtConnector(
                cache_dir=self.config.cache_dir / "uniprot",
                rate_limit=self.config.uniprot.rate_limit,
            )

        if self.config.pubchem.enabled:
            self.connectors["pubchem"] = PubChemConnector(
                cache_dir=self.config.cache_dir / "pubchem",
                rate_limit=self.config.pubchem.rate_limit,
            )

        if self.config.zinc.enabled:
            self.connectors["zinc"] = ZINCConnector(
                cache_dir=self.config.cache_dir / "zinc",
            )

        if self.config.uspto.enabled:
            self.connectors["uspto"] = USPTOConnector(
                cache_dir=self.config.cache_dir / "uspto",
            )

        if self.config.retrobiocat.enabled and self.config.retrobiocat.path:
            self.connectors["retrobiocat"] = RetroBioCatConnector(
                rxns_yaml_path=self.config.retrobiocat.path,
                cache_dir=self.config.cache_dir / "retrobiocat",
            )

        # Connect all
        for name, connector in self.connectors.items():
            try:
                if connector.connect():
                    logger.info(f"✓ Connected to {name}")
                else:
                    logger.warning(f"✗ Failed to connect to {name}")
            except Exception as e:
                logger.error(f"✗ Error connecting to {name}: {e}")

    def disconnect_all(self) -> None:
        """Disconnect all connectors."""
        for name, connector in self.connectors.items():
            try:
                connector.disconnect()
                logger.info(f"Disconnected from {name}")
            except Exception as e:
                logger.error(f"Error disconnecting from {name}: {e}")

    def build_from_local_sources(self) -> Iterator[Dict[str, Any]]:
        """Build knowledge base from local data sources (BKMS, BRENDA, EnzyExtract)."""
        logger.info("Processing local data sources...")

        # Process in order: BKMS as core, then enrich with kinetics
        for source_name in ["bkms", "brenda", "enzyextract"]:
            if source_name not in self.connectors:
                continue

            connector = self.connectors[source_name]
            if not connector.is_available():
                logger.warning(f"Skipping {source_name} (not available)")
                continue

            logger.info(f"Processing {source_name}...")

            try:
                for record in tqdm(
                    connector.get_all_records(),
                    desc=f"Reading {source_name}",
                ):
                    self.stats["total_records"] += 1
                    self.stats["by_source"][source_name] += 1

                    # Track EC numbers
                    for ec in record.get("ec_numbers", []):
                        if self.validator.validate_ec_number(ec):
                            self.stats["by_ec"][ec] += 1

                    yield record

            except Exception as e:
                error_msg = f"Error processing {source_name}: {e}"
                logger.error(error_msg)
                self.stats["errors"].append(error_msg)

    def enrich_with_api_sources(
        self,
        records: List[Dict[str, Any]],
        sources: Optional[List[str]] = None,
    ) -> Iterator[Dict[str, Any]]:
        """Enrich records with data from API sources."""
        sources = sources or ["kegg", "uniprot", "pubchem"]
        logger.info(f"Enriching with API sources: {sources}")

        for record in tqdm(records, desc="Enriching records"):
            enriched = record.copy()

            # Get EC number for API queries
            primary_ec = record.get("primary_ec")
            if not primary_ec:
                yield enriched
                continue

            # Query each API source
            for source_name in sources:
                if source_name not in self.connectors:
                    continue

                connector = self.connectors[source_name]
                if not connector.is_available():
                    continue

                try:
                    result = connector.query_by_ec(primary_ec)
                    if result.success and result.data:
                        # Merge data
                        enriched = self._merge_records(enriched, result.data, source_name)
                except Exception as e:
                    logger.warning(f"Failed to query {source_name} for EC {primary_ec}: {e}")

            yield enriched

    def _merge_records(
        self,
        base: Dict[str, Any],
        api_data: Any,
        source: str,
    ) -> Dict[str, Any]:
        """Merge API data into base record."""
        merged = base.copy()

        if not api_data:
            return merged

        # Handle different data types
        if isinstance(api_data, list):
            api_record = api_data[0] if api_data else {}
        else:
            api_record = api_data

        # Merge specific fields based on source
        if source == "kegg":
            if "pathway" not in merged:
                merged["pathway"] = {}
            merged["pathway"].update({
                "kegg_pathways": api_record.get("kegg_specific", {}).get("linked_pathways", []),
                "kegg_reactions": api_record.get("kegg_specific", {}).get("linked_reactions", []),
            })

        elif source == "uniprot":
            # Add UniProt IDs if not present
            uniprot_ids = api_record.get("enzyme", {}).get("uniprot_ids", [])
            if uniprot_ids and not merged.get("enzyme", {}).get("uniprot_ids"):
                if "enzyme" not in merged:
                    merged["enzyme"] = {}
                merged["enzyme"]["uniprot_ids"] = uniprot_ids

            # Add sequence if available
            sequence = api_record.get("enzyme", {}).get("sequence")
            if sequence and not merged.get("enzyme", {}).get("sequence"):
                merged["enzyme"]["sequence"] = sequence

        elif source == "pubchem":
            # Add PubChem compound data
            substrates = merged.get("reaction", {}).get("substrates", [])
            for substrate in substrates:
                if not substrate.get("cid"):
                    compound_data = api_record.get("compound", {})
                    if compound_data:
                        substrate["cid"] = str(compound_data.get("cid", ""))
                        substrate["smiles"] = compound_data.get("smiles", "")

        return merged

    def build_ec_index(
        self,
        records: Iterator[Dict[str, Any]],
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Build EC number index from records."""
        logger.info("Building EC number index...")
        ec_index: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

        for record in tqdm(records, desc="Indexing by EC"):
            for ec in record.get("ec_numbers", []):
                if self.validator.validate_ec_number(ec):
                    ec_index[ec].append(record)

        logger.info(f"Indexed {len(ec_index)} unique EC numbers")
        return dict(ec_index)

    def save_to_jsonl(
        self,
        records: Iterator[Dict[str, Any]],
        output_path: Path,
        compress: bool = False,
    ) -> None:
        """Save records to JSONL file."""
        logger.info(f"Saving to {output_path}...")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        if compress:
            import gzip
            opener = gzip.open
            output_path = output_path.with_suffix(output_path.suffix + ".gz")
        else:
            opener = open

        count = 0
        with opener(output_path, "wt", encoding="utf-8") as f:
            for record in tqdm(records, desc="Writing records"):
                f.write(json.dumps(record, ensure_ascii=False) + "\n")
                count += 1

        logger.info(f"✓ Saved {count} records to {output_path}")

    def save_to_parquet(
        self,
        records: Iterator[Dict[str, Any]],
        output_path: Path,
    ) -> None:
        """Save records to Parquet file."""
        logger.info(f"Saving to {output_path}...")

        import pandas as pd

        # Collect all records
        record_list = list(tqdm(records, desc="Collecting records"))

        if not record_list:
            logger.warning("No records to save")
            return

        # Convert to DataFrame
        df = pd.DataFrame(record_list)

        # Save
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(output_path, index=False, compression="snappy")

        logger.info(f"✓ Saved {len(df)} records to {output_path}")

    def build_full_knowledge_base(self, enrich: bool = False) -> Path:
        """Build complete knowledge base from all sources."""
        logger.info("=" * 60)
        logger.info("Building unified knowledge base")
        logger.info("=" * 60)

        self.initialize_connectors()

        # Collect records from local sources
        records = list(self.build_from_local_sources())

        # Optionally enrich with API data
        if enrich:
            records = list(self.enrich_with_api_sources(records))

        # Save to output
        output_format = self.config.output_format
        output_file = self.config.output_dir / f"knowledge_base.{output_format}"

        if output_format == "jsonl":
            self.save_to_jsonl(
                iter(records),
                output_file,
                compress=self.config.compress_output,
            )
        elif output_format == "parquet":
            self.save_to_parquet(iter(records), output_file)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")

        # Save stats
        stats_file = self.config.output_dir / "build_stats.json"
        with open(stats_file, "w") as f:
            json.dump(self.stats, f, indent=2)

        logger.info("=" * 60)
        logger.info(f"Build complete! Output: {output_file}")
        logger.info(f"Total records: {self.stats['total_records']}")
        logger.info(f"Unique EC numbers: {len(self.stats['by_ec'])}")
        logger.info(f"Stats saved to: {stats_file}")
        logger.info("=" * 60)

        self.disconnect_all()

        return output_file

    def __enter__(self) -> "KnowledgeBaseBuilder":
        self.initialize_connectors()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.disconnect_all()
