from __future__ import annotations

import ast
import math
from typing import Any, Dict, Iterable, List, Optional

Number = float | int


def _is_missing(value: Any) -> bool:
    """Return True for None/NaN/empty strings."""
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str) and not value.strip():
        return True
    return False


def _ensure_list(value: Any) -> List[Any]:
    """Best-effort conversion of semi-structured columns into a list."""
    if _is_missing(value):
        return []
    if isinstance(value, list):
        return [v for v in value if not _is_missing(v)]
    if isinstance(value, tuple):
        return [v for v in value if not _is_missing(v)]
    if hasattr(value, "tolist"):
        try:
            as_list = value.tolist()
            if isinstance(as_list, list):
                return [v for v in as_list if not _is_missing(v)]
        except Exception:
            pass
    if isinstance(value, str):
        text = value.strip()
        # Try to parse simple repr lists like "['1.1.1.1', '2.2.2.2']"
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = ast.literal_eval(text)
                if isinstance(parsed, list):
                    return [v for v in parsed if not _is_missing(v)]
            except (ValueError, SyntaxError):
                pass
        # Fallback: split on commas/semicolons
        parts = [p.strip() for p in text.replace(";", ",").split(",")]
        return [p for p in parts if p]
    return [value]


def _to_float(value: Any) -> Optional[float]:
    if _is_missing(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _standardize_kcat_unit(raw_text: Any, default: str = "s^-1") -> Optional[str]:
    if _is_missing(raw_text):
        return default
    text = str(raw_text).lower()
    if any(token in text for token in ["s^-1", "s-1", "/s", "per second"]):
        return "s^-1"
    return default


def _standardize_km_unit(raw_text: Any, default: str = "mM") -> Optional[str]:
    if _is_missing(raw_text):
        return default
    text = str(raw_text).lower()
    if "µm" in text or "um" in text:
        return "uM"
    if "mm" in text:
        return "mM"
    if "m" in text:
        return "M"
    return default


def _parse_cofactors(value: Any) -> List[Dict[str, Any]]:
    cofactors: List[Dict[str, Any]] = []
    for name in _ensure_list(value):
        cofactors.append(
            {
                "name": str(name),
                "role": "cofactor",
                "stoichiometry": 1.0,
                "smiles": None,
            }
        )
    return cofactors


def normalize_list(value: Any) -> List[Any]:
    """Public wrapper for list normalization."""
    return _ensure_list(value)


def parse_cofactors(value: Any) -> List[Dict[str, Any]]:
    """Public wrapper to translate cofactor strings into structured entries."""
    return _parse_cofactors(value)


def build_base_record(
    dataset: str, row_index: int | str, raw_ids: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Skeleton record with required top-level keys and defaults."""
    rid = f"{dataset}:{row_index}"
    return {
        "id": rid,
        "source": {
            "dataset": dataset,
            "row_index": row_index,
            "raw_ids": {k: v for k, v in (raw_ids or {}).items() if not _is_missing(v)},
        },
        "ec_numbers": [],
        "primary_ec": None,
        "enzyme": {
            "name": None,
            "name_full": None,
            "uniprot_ids": [],
            "sequence": None,
            "pdb_ids": [],
            "organism": None,
            "mutant": None,
            "mutation_flag": False,
        },
        "reaction": {
            "direction": "unknown",
            "equation_text": None,
            "substrates": [],
            "products": [],
            "cofactors": [],
            "reaction_smarts": None,
            "mechanism_type": None,
        },
        "kinetics": {
            "kcat": {
                "value": None,
                "unit": None,
                "raw_text": None,
                "condition_id": None,
            },
            "km": {
                "value": None,
                "unit": None,
                "raw_text": None,
                "condition_id": None,
            },
            "kcat_over_km": {
                "value": None,
                "unit": None,
            },
        },
        "conditions": {
            "id": None,
            "pH": None,
            "temperature": None,
            "temperature_unit": None,
            "ionic_strength": None,
            "buffer": None,
            "comments": None,
        },
        "text": {
            "descriptor": None,
            "canonical_description": None,
            "notes": None,
        },
    }


def parse_bkms_equation(equation: str) -> Dict[str, Any]:
    """
    Split a BKMS reaction equation into substrates/products while preserving order.
    Supports separators: '<=>', '=>', '<=', '->', '='.
    """
    if not equation:
        return {"direction": "unknown", "substrates": [], "products": [], "equation_text": None}

    direction = "unknown"
    separators = ["<=>", "=>", "<=", "->", "="]
    left, right = equation, ""
    for sep in separators:
        if sep in equation:
            left, right = equation.split(sep, 1)
            if sep in ("=>", "->", "="):
                direction = "forward"
            elif sep == "<=":
                direction = "reverse"
            else:
                direction = "unknown"
            break

    def split_side(side: str, role: str) -> List[Dict[str, Any]]:
        parts = [p.strip() for p in side.split("+") if p.strip()]
        result: List[Dict[str, Any]] = []
        for part in parts:
            stoich = 1.0
            # Capture stoichiometry prefix like "2 H2O"
            tokens = part.split()
            if tokens and tokens[0].replace(".", "", 1).isdigit():
                try:
                    stoich = float(tokens[0])
                    part = " ".join(tokens[1:]) or tokens[0]
                except ValueError:
                    pass
            result.append(
                {
                    "name": part,
                    "name_full": None,
                    "role": role,
                    "stoichiometry": stoich,
                    "smiles": None,
                    "cid": None,
                    "brenda_id": None,
                }
            )
        return result

    substrates = split_side(left, "substrate")
    products = split_side(right, "product")

    return {
        "direction": direction,
        "substrates": substrates,
        "products": products,
        "equation_text": equation,
    }


def build_conditions(
    dataset: str,
    row_index: int | str,
    pH: Any = None,
    temperature: Any = None,
    temperature_unit: Optional[str] = "C",
    ionic_strength: Any = None,
    buffer: Any = None,
    comments: Any = None,
) -> Dict[str, Any]:
    condition_id = f"cond_{dataset}_{row_index}"
    return {
        "id": condition_id,
        "pH": _to_float(pH),
        "temperature": _to_float(temperature),
        "temperature_unit": temperature_unit if not _is_missing(temperature) else None,
        "ionic_strength": ionic_strength if not _is_missing(ionic_strength) else None,
        "buffer": buffer if not _is_missing(buffer) else None,
        "comments": comments if not _is_missing(comments) else None,
    }


def populate_enzyme_section(
    record: Dict[str, Any],
    *,
    name: Any = None,
    name_full: Any = None,
    organism: Any = None,
    mutant: Any = None,
    uniprot_ids: Iterable[Any] | None = None,
    sequence: Any = None,
    pdb_ids: Iterable[Any] | None = None,
) -> None:
    record["enzyme"]["name"] = None if _is_missing(name) else name
    record["enzyme"]["name_full"] = None if _is_missing(name_full) else name_full
    record["enzyme"]["organism"] = None if _is_missing(organism) else organism
    record["enzyme"]["mutant"] = None if _is_missing(mutant) else mutant
    record["enzyme"]["mutation_flag"] = not _is_missing(mutant) and str(mutant).upper() != "WT"
    record["enzyme"]["uniprot_ids"] = [
        uid for uid in (uniprot_ids or []) if not _is_missing(uid)
    ]
    record["enzyme"]["sequence"] = None if _is_missing(sequence) else sequence
    record["enzyme"]["pdb_ids"] = [pid for pid in (pdb_ids or []) if not _is_missing(pid)]


def populate_reaction(
    record: Dict[str, Any],
    *,
    direction: str = "unknown",
    equation_text: Any = None,
    substrates: List[Dict[str, Any]] | None = None,
    products: List[Dict[str, Any]] | None = None,
    cofactors: List[Dict[str, Any]] | None = None,
) -> None:
    record["reaction"]["direction"] = direction
    record["reaction"]["equation_text"] = None if _is_missing(equation_text) else equation_text
    record["reaction"]["substrates"] = substrates or []
    record["reaction"]["products"] = products or []
    record["reaction"]["cofactors"] = cofactors or []


def populate_kinetics(
    record: Dict[str, Any],
    *,
    dataset: str,
    row_index: int | str,
    kcat_value: Any = None,
    kcat_raw: Any = None,
    km_value: Any = None,
    km_raw: Any = None,
    kcat_km_value: Any = None,
) -> None:
    cond_id = f"cond_{dataset}_{row_index}"
    record["kinetics"]["kcat"] = {
        "value": _to_float(kcat_value),
        "unit": _standardize_kcat_unit(kcat_raw),
        "raw_text": None if _is_missing(kcat_raw) else str(kcat_raw),
        "condition_id": cond_id,
    }
    record["kinetics"]["km"] = {
        "value": _to_float(km_value),
        "unit": _standardize_km_unit(km_raw),
        "raw_text": None if _is_missing(km_raw) else str(km_raw),
        "condition_id": cond_id,
    }
    record["kinetics"]["kcat_over_km"] = {
        "value": _to_float(kcat_km_value),
        "unit": "s^-1 mM^-1" if not _is_missing(kcat_km_value) else None,
    }


def populate_text(
    record: Dict[str, Any],
    *,
    descriptor: Any = None,
    canonical_description: Any = None,
    notes: Any = None,
) -> None:
    record["text"]["descriptor"] = None if _is_missing(descriptor) else descriptor
    record["text"]["canonical_description"] = (
        None if _is_missing(canonical_description) else canonical_description
    )
    record["text"]["notes"] = None if _is_missing(notes) else notes


__all__ = [
    "build_base_record",
    "build_conditions",
    "normalize_list",
    "parse_cofactors",
    "parse_bkms_equation",
    "populate_enzyme_section",
    "populate_reaction",
    "populate_kinetics",
    "populate_text",
]
