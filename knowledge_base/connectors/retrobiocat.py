"""RetroBioCat biocatalysis database connector."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class RetroBioCatConnector(BaseConnector):
    """Connector for RetroBioCat reaction templates and activity data."""

    SOURCE_NAME = "RetroBioCat"
    ENTITY_TYPE = "reaction_template"

    def __init__(
        self,
        rxns_yaml_path: Optional[Path] = None,
        activity_data_path: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        **kwargs,
    ):
        super().__init__(cache_dir=cache_dir, **kwargs)
        self.rxns_yaml_path = rxns_yaml_path
        self.activity_data_path = activity_data_path
        self._templates: List[Dict[str, Any]] = []

    def connect(self) -> bool:
        """Load RetroBioCat data files."""
        if self.rxns_yaml_path and self.rxns_yaml_path.exists():
            try:
                import yaml
                with open(self.rxns_yaml_path, "r") as f:
                    data = yaml.safe_load(f)
                    if isinstance(data, dict):
                        self._templates = list(data.values())
                    logger.info(f"Loaded {len(self._templates)} RetroBioCat templates")
                return True
            except Exception as e:
                logger.error(f"Failed to load RetroBioCat YAML: {e}")
                return False
        logger.warning("RetroBioCat data path not configured")
        return False

    def disconnect(self) -> None:
        """Clear loaded data."""
        self._templates = []

    def is_available(self) -> bool:
        """Check if data is loaded."""
        return len(self._templates) > 0

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query templates by EC number (not directly supported)."""
        return ConnectorResult(
            success=True,
            data=[],
            source=self.SOURCE_NAME,
            metadata={"note": "EC queries not directly supported"},
        )

    def query_by_compound(self, identifier: str, id_type: str = "name") -> ConnectorResult:
        """Query templates involving compound."""
        records = []
        identifier_lower = identifier.lower()

        for template in self._templates:
            # Check if compound appears in template
            rxn_str = str(template.get("reaction_smarts", "")).lower()
            if identifier_lower in rxn_str:
                records.append(self.to_unified_schema(template))

        return ConnectorResult(
            success=True,
            data=records,
            source=self.SOURCE_NAME,
            metadata={"count": len(records)},
        )

    def query_by_reaction(
        self, substrates: Optional[List[str]] = None, products: Optional[List[str]] = None
    ) -> ConnectorResult:
        """Query templates by reaction components."""
        return ConnectorResult(success=True, data=[], source=self.SOURCE_NAME)

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over all templates."""
        for idx, template in enumerate(self._templates):
            yield self.to_unified_schema(template, idx)

    def to_unified_schema(self, record: Dict[str, Any], idx: Optional[int] = None) -> Dict[str, Any]:
        """Convert RetroBioCat template to unified schema."""
        idx = idx or 0

        return {
            "id": f"RetroBioCat:{idx}",
            "source": {
                "dataset": "RetroBioCat",
                "row_index": idx,
            },
            "ec_numbers": [],
            "primary_ec": None,
            "enzyme": {
                "name": record.get("enzyme_type"),
                "enzyme_class": record.get("type"),
            },
            "reaction": {
                "reaction_smarts": record.get("reaction_smarts"),
                "template_name": record.get("name"),
                "cofactors": record.get("cofactors", []),
            },
            "retrosynthesis": {
                "template": record.get("reaction_smarts"),
                "enzyme_type": record.get("enzyme_type"),
                "positive_tests": record.get("positive_tests", 0),
                "negative_tests": record.get("negative_tests", 0),
            },
        }
