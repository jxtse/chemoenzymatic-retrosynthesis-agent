"""ZINC purchasable compounds database connector."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class ZINCConnector(BaseConnector):
    """Connector for ZINC database (stub implementation)."""

    SOURCE_NAME = "ZINC"
    ENTITY_TYPE = "compound"

    def __init__(self, cache_dir: Optional[Path] = None, **kwargs):
        super().__init__(cache_dir=cache_dir, **kwargs)

    def connect(self) -> bool:
        """Connect to ZINC (placeholder)."""
        logger.warning("ZINC connector is a stub implementation")
        return True

    def disconnect(self) -> None:
        """Disconnect."""
        pass

    def is_available(self) -> bool:
        """Check availability."""
        return False

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Not supported."""
        return ConnectorResult(
            success=False,
            error="ZINC does not support EC queries",
            source=self.SOURCE_NAME,
        )

    def query_by_compound(self, identifier: str, id_type: str = "name") -> ConnectorResult:
        """Stub implementation."""
        return ConnectorResult(
            success=False,
            error="ZINC connector not fully implemented",
            source=self.SOURCE_NAME,
        )

    def query_by_reaction(
        self, substrates: Optional[List[str]] = None, products: Optional[List[str]] = None
    ) -> ConnectorResult:
        """Not supported."""
        return ConnectorResult(success=False, error="Not implemented", source=self.SOURCE_NAME)

    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Not supported."""
        return iter([])

    def to_unified_schema(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """Convert to unified schema."""
        return record
