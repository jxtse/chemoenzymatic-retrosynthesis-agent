"""Base connector interface for all data sources."""

from __future__ import annotations

import hashlib
import json
import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")


@dataclass
class ConnectorResult:
    """Standardized result from a connector query."""

    success: bool
    data: Optional[Any] = None
    error: Optional[str] = None
    source: str = ""
    query_time: float = 0.0
    cached: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CacheEntry:
    """Cache entry with timestamp."""

    data: Any
    timestamp: datetime
    ttl_seconds: int = 3600


class BaseConnector(ABC):
    """Abstract base class for all database connectors."""

    SOURCE_NAME: str = "base"
    ENTITY_TYPE: str = "record"  # enzyme, compound, reaction, etc.

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        rate_limit: float = 1.0,
        batch_size: int = 100,
    ):
        self.cache_dir = cache_dir
        self.rate_limit = rate_limit
        self.batch_size = batch_size
        self._last_request_time: float = 0
        self._memory_cache: Dict[str, CacheEntry] = {}

        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

    def _rate_limit_wait(self) -> None:
        """Enforce rate limiting between requests."""
        if self.rate_limit > 0:
            elapsed = time.time() - self._last_request_time
            min_interval = 1.0 / self.rate_limit
            if elapsed < min_interval:
                time.sleep(min_interval - elapsed)
        self._last_request_time = time.time()

    def _cache_key(self, query: str, params: Optional[Dict[str, Any]] = None) -> str:
        """Generate cache key from query and params."""
        key_str = f"{self.SOURCE_NAME}:{query}"
        if params:
            key_str += ":" + json.dumps(params, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()

    def _get_from_cache(self, key: str) -> Optional[Any]:
        """Try to get result from cache."""
        # Check memory cache first
        if key in self._memory_cache:
            entry = self._memory_cache[key]
            age = (datetime.now() - entry.timestamp).total_seconds()
            if age < entry.ttl_seconds:
                return entry.data

        # Check disk cache
        if self.cache_dir:
            cache_file = self.cache_dir / f"{key}.json"
            if cache_file.exists():
                try:
                    with open(cache_file, "r", encoding="utf-8") as f:
                        cached = json.load(f)
                    if time.time() - cached.get("timestamp", 0) < cached.get("ttl", 3600):
                        return cached.get("data")
                except (json.JSONDecodeError, IOError):
                    pass
        return None

    def _save_to_cache(self, key: str, data: Any, ttl_seconds: int = 3600) -> None:
        """Save result to cache."""
        # Save to memory cache
        self._memory_cache[key] = CacheEntry(
            data=data,
            timestamp=datetime.now(),
            ttl_seconds=ttl_seconds,
        )

        # Save to disk cache
        if self.cache_dir:
            cache_file = self.cache_dir / f"{key}.json"
            try:
                with open(cache_file, "w", encoding="utf-8") as f:
                    json.dump({
                        "data": data,
                        "timestamp": time.time(),
                        "ttl": ttl_seconds,
                    }, f)
            except (TypeError, IOError) as e:
                logger.warning(f"Failed to save to disk cache: {e}")

    @abstractmethod
    def connect(self) -> bool:
        """Establish connection to the data source."""
        pass

    @abstractmethod
    def disconnect(self) -> None:
        """Close connection to the data source."""
        pass

    @abstractmethod
    def is_available(self) -> bool:
        """Check if the data source is available."""
        pass

    @abstractmethod
    def query_by_ec(self, ec_number: str) -> ConnectorResult:
        """Query records by EC number."""
        pass

    @abstractmethod
    def query_by_compound(
        self,
        identifier: str,
        id_type: str = "name"  # name, smiles, inchi, cid, etc.
    ) -> ConnectorResult:
        """Query records by compound identifier."""
        pass

    @abstractmethod
    def query_by_reaction(
        self,
        substrates: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> ConnectorResult:
        """Query records by reaction components."""
        pass

    @abstractmethod
    def get_all_records(self) -> Iterator[Dict[str, Any]]:
        """Iterate over all records in the data source."""
        pass

    @abstractmethod
    def to_unified_schema(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """Convert a record to the unified schema format."""
        pass

    def batch_query_by_ec(self, ec_numbers: List[str]) -> List[ConnectorResult]:
        """Query multiple EC numbers in batch."""
        results = []
        for ec in ec_numbers:
            results.append(self.query_by_ec(ec))
        return results

    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about the data source."""
        return {
            "source": self.SOURCE_NAME,
            "entity_type": self.ENTITY_TYPE,
            "available": self.is_available(),
            "cache_entries": len(self._memory_cache),
        }

    def __enter__(self) -> "BaseConnector":
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.disconnect()
