"""Unified schema definition for the knowledge base."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class UnifiedRecord:
    """Unified record structure for all database sources."""

    id: str
    source: Dict[str, Any]
    ec_numbers: List[str] = field(default_factory=list)
    primary_ec: Optional[str] = None

    enzyme: Dict[str, Any] = field(default_factory=dict)
    reaction: Dict[str, Any] = field(default_factory=dict)
    kinetics: Dict[str, Any] = field(default_factory=dict)
    conditions: Dict[str, Any] = field(default_factory=dict)
    text: Dict[str, Any] = field(default_factory=dict)

    # Optional fields for specific data types
    compound: Optional[Dict[str, Any]] = None
    pathway: Optional[Dict[str, Any]] = None
    protein: Optional[Dict[str, Any]] = None
    patent: Optional[Dict[str, Any]] = None

    # Source-specific extra fields
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "id": self.id,
            "source": self.source,
            "ec_numbers": self.ec_numbers,
            "primary_ec": self.primary_ec,
            "enzyme": self.enzyme,
            "reaction": self.reaction,
            "kinetics": self.kinetics,
            "conditions": self.conditions,
            "text": self.text,
        }

        if self.compound is not None:
            result["compound"] = self.compound
        if self.pathway is not None:
            result["pathway"] = self.pathway
        if self.protein is not None:
            result["protein"] = self.protein
        if self.patent is not None:
            result["patent"] = self.patent
        if self.extra:
            result["extra"] = self.extra

        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "UnifiedRecord":
        """Create from dictionary."""
        return cls(
            id=data.get("id", ""),
            source=data.get("source", {}),
            ec_numbers=data.get("ec_numbers", []),
            primary_ec=data.get("primary_ec"),
            enzyme=data.get("enzyme", {}),
            reaction=data.get("reaction", {}),
            kinetics=data.get("kinetics", {}),
            conditions=data.get("conditions", {}),
            text=data.get("text", {}),
            compound=data.get("compound"),
            pathway=data.get("pathway"),
            protein=data.get("protein"),
            patent=data.get("patent"),
            extra=data.get("extra", {}),
        )


class SchemaValidator:
    """Validator for unified schema."""

    @staticmethod
    def validate_ec_number(ec: str) -> bool:
        """Validate EC number format (e.g., 1.2.3.4)."""
        if not ec:
            return False
        parts = ec.split(".")
        if len(parts) != 4:
            return False
        for part in parts:
            if part == "-":
                continue
            try:
                int(part)
            except ValueError:
                return False
        return True

    @staticmethod
    def validate_smiles(smiles: str) -> bool:
        """Validate SMILES string (basic check)."""
        if not smiles:
            return False
        # Basic check - real validation requires RDKit
        invalid_chars = set("[]()<>=# ") - set(smiles)
        return len(smiles) > 0

    @staticmethod
    def normalize_kinetic_value(value: Any, unit: str) -> Optional[float]:
        """Normalize kinetic values to standard units."""
        if value is None:
            return None

        try:
            val = float(value)
        except (TypeError, ValueError):
            return None

        # Normalize based on unit
        unit_lower = unit.lower() if unit else ""

        # kcat to s^-1
        if "min" in unit_lower and "kcat" in unit_lower:
            return val / 60.0

        # Km to mM
        if "um" in unit_lower or "µm" in unit_lower:
            return val / 1000.0
        elif unit_lower == "m":
            return val * 1000.0

        return val
