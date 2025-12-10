"""Configuration management for the unified knowledge base."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


@dataclass
class DataSourceConfig:
    """Configuration for a single data source."""

    enabled: bool = True
    path: Optional[Path] = None
    api_url: Optional[str] = None
    api_key: Optional[str] = None
    cache_dir: Optional[Path] = None
    batch_size: int = 1000
    rate_limit: float = 1.0  # requests per second
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class KnowledgeBaseConfig:
    """Configuration for the unified knowledge base."""

    # Base paths
    data_dir: Path = field(default_factory=lambda: Path("data"))
    output_dir: Path = field(default_factory=lambda: Path("knowledge_base_output"))
    cache_dir: Path = field(default_factory=lambda: Path(".cache"))

    # Individual data source configs
    bkms: DataSourceConfig = field(default_factory=DataSourceConfig)
    brenda: DataSourceConfig = field(default_factory=DataSourceConfig)
    enzyextract: DataSourceConfig = field(default_factory=DataSourceConfig)
    kegg: DataSourceConfig = field(default_factory=DataSourceConfig)
    uniprot: DataSourceConfig = field(default_factory=DataSourceConfig)
    pubchem: DataSourceConfig = field(default_factory=DataSourceConfig)
    zinc: DataSourceConfig = field(default_factory=DataSourceConfig)
    uspto: DataSourceConfig = field(default_factory=DataSourceConfig)
    retrobiocat: DataSourceConfig = field(default_factory=DataSourceConfig)

    # Processing options
    parallel_workers: int = 4
    chunk_size: int = 10000
    normalize_smiles: bool = True
    validate_ec_numbers: bool = True

    # Output options
    output_format: str = "jsonl"  # jsonl, parquet, sqlite
    compress_output: bool = False

    @classmethod
    def from_yaml(cls, path: Path) -> "KnowledgeBaseConfig":
        """Load configuration from YAML file."""
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "KnowledgeBaseConfig":
        """Create config from dictionary."""
        config = cls()

        # Set base paths
        if "data_dir" in data:
            config.data_dir = Path(data["data_dir"])
        if "output_dir" in data:
            config.output_dir = Path(data["output_dir"])
        if "cache_dir" in data:
            config.cache_dir = Path(data["cache_dir"])

        # Set data source configs
        source_names = [
            "bkms", "brenda", "enzyextract", "kegg",
            "uniprot", "pubchem", "zinc", "uspto", "retrobiocat"
        ]
        for name in source_names:
            if name in data:
                src_data = data[name]
                # Handle relative paths by prepending data_dir
                src_path = None
                if src_data.get("path"):
                    src_path = Path(src_data["path"])
                    # If path is relative and doesn't exist, try prepending data_dir
                    if not src_path.is_absolute() and not src_path.exists():
                        src_path = config.data_dir / src_data["path"]
                src_config = DataSourceConfig(
                    enabled=src_data.get("enabled", True),
                    path=src_path,
                    api_url=src_data.get("api_url"),
                    api_key=src_data.get("api_key") or os.environ.get(f"{name.upper()}_API_KEY"),
                    cache_dir=Path(src_data["cache_dir"]) if src_data.get("cache_dir") else None,
                    batch_size=src_data.get("batch_size", 1000),
                    rate_limit=src_data.get("rate_limit", 1.0),
                    extra=src_data.get("extra", {}),
                )
                setattr(config, name, src_config)

        # Set processing options
        if "parallel_workers" in data:
            config.parallel_workers = data["parallel_workers"]
        if "chunk_size" in data:
            config.chunk_size = data["chunk_size"]
        if "normalize_smiles" in data:
            config.normalize_smiles = data["normalize_smiles"]
        if "validate_ec_numbers" in data:
            config.validate_ec_numbers = data["validate_ec_numbers"]

        # Set output options
        if "output_format" in data:
            config.output_format = data["output_format"]
        if "compress_output" in data:
            config.compress_output = data["compress_output"]

        return config

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        def source_to_dict(src: DataSourceConfig) -> Dict[str, Any]:
            return {
                "enabled": src.enabled,
                "path": str(src.path) if src.path else None,
                "api_url": src.api_url,
                "api_key": "***" if src.api_key else None,  # Don't expose API keys
                "cache_dir": str(src.cache_dir) if src.cache_dir else None,
                "batch_size": src.batch_size,
                "rate_limit": src.rate_limit,
                "extra": src.extra,
            }

        return {
            "data_dir": str(self.data_dir),
            "output_dir": str(self.output_dir),
            "cache_dir": str(self.cache_dir),
            "bkms": source_to_dict(self.bkms),
            "brenda": source_to_dict(self.brenda),
            "enzyextract": source_to_dict(self.enzyextract),
            "kegg": source_to_dict(self.kegg),
            "uniprot": source_to_dict(self.uniprot),
            "pubchem": source_to_dict(self.pubchem),
            "zinc": source_to_dict(self.zinc),
            "uspto": source_to_dict(self.uspto),
            "retrobiocat": source_to_dict(self.retrobiocat),
            "parallel_workers": self.parallel_workers,
            "chunk_size": self.chunk_size,
            "normalize_smiles": self.normalize_smiles,
            "validate_ec_numbers": self.validate_ec_numbers,
            "output_format": self.output_format,
            "compress_output": self.compress_output,
        }

    def save_yaml(self, path: Path) -> None:
        """Save configuration to YAML file."""
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    @classmethod
    def default(cls) -> "KnowledgeBaseConfig":
        """Create default configuration with common paths."""
        config = cls()

        # Set default paths for local data files
        base = Path(".")
        config.bkms.path = base / "Reactions_BKMS.csv"
        config.brenda.path = base / "brenda_kcat_v3.parquet"
        config.enzyextract.path = base / "EnzyExtractDB_176463.parquet"

        # Set API URLs
        config.kegg.api_url = "https://rest.kegg.jp"
        config.uniprot.api_url = "https://rest.uniprot.org"
        config.pubchem.api_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        config.zinc.api_url = "https://zinc.docking.org/api"
        config.brenda.api_url = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

        return config
