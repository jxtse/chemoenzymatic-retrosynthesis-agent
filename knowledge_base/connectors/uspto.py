"""USPTO patent database connector."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from .base import BaseConnector, ConnectorResult

logger = logging.getLogger(__name__)


class USPTOConnector(BaseConnector):
    """Connector for USPTO patent database (stub implementation)."""

    SOURCE_NAME = "USPTO"
    ENTITY_TYPE = "patent"

    def __init__(self, cache_dir: Optional[Path] = None, **kwargs):
        super().__init__(cache_dir=cache_dir, **kwargs)

    def connect(self) -> bool:
        """Connect to USPTO (placeholder)."""
        logger.warning("USPTO connector is a stub implementation")
        return True

    def disconnect(self) -> None:
        """Disconnect."""
        pass

    def is_available(self) -> bool:
        """Check availability."""
        return False

    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Search patents by EC number."""
        return ConnectorResult(
            success=False,
            error="USPTO connector not fully implemented",
            source=self.SOURCE_NAME,
        )

    def query_by_compound(self, identifier: str, id_type: str = "name") -> ConnectorResult:
        """Search patents by compound."""
        return ConnectorResult(
            success=False,
            error="USPTO connector not fully implemented",
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
