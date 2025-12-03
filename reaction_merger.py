from __future__ import annotations

import argparse
import copy
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Optional

from reaction_unifier import bkms_records, brenda_records, enzyextract_records


def _index_by_ec(records: Iterable[Dict[str, object]]) -> DefaultDict[str, List[Dict[str, object]]]:
    """Build a mapping from EC number to list of measurement records."""
    index: DefaultDict[str, List[Dict[str, object]]] = defaultdict(list)
    for rec in records:
        for ec in rec.get("ec_numbers", []) or []:
            index[str(ec)].append(rec)
    return index


def merge_reactions(
    bkms_path: Path,
    enzy_path: Path | None,
    brenda_path: Path | None,
    limit: int | None = None,
    view: str = "full",
) -> Iterator[Dict[str, object]]:
    """Yield BKMS reaction cores enriched with matching kinetics by EC number."""
    enzy_index: DefaultDict[str, List[Dict[str, object]]] = defaultdict(list)
    brenda_index: DefaultDict[str, List[Dict[str, object]]] = defaultdict(list)

    if enzy_path:
        enzy_index = _index_by_ec(enzyextract_records(enzy_path))
    if brenda_path:
        brenda_index = _index_by_ec(brenda_records(brenda_path))

    count = 0
    for core in bkms_records(bkms_path):
        if limit is not None and count >= limit:
            break
        count += 1

        ecs = core.get("ec_numbers", []) or []
        matched: List[Dict[str, object]] = []
        for ec in ecs:
            matched.extend(enzy_index.get(ec, []))
            matched.extend(brenda_index.get(ec, []))

        # Remove duplicates by id while preserving order
        seen_ids = set()
        uniq: List[Dict[str, object]] = []
        for rec in matched:
            rec_id = rec.get("id")
            if rec_id in seen_ids:
                continue
            seen_ids.add(rec_id)
            uniq.append(annotate_measurement(rec, core.get("reaction", {}) or {}))

        enriched = dict(core)
        enriched["measurements"] = uniq
        if view == "minimal":
            yield simplify_record(enriched)
        elif view == "flat":
            for row in flatten_record(enriched):
                yield row
        else:
            yield enriched


def _names_from_compounds(items: List[Dict[str, object]]) -> List[str]:
    names: List[str] = []
    for itm in items:
        name = itm.get("name")
        if name:
            names.append(str(name))
    return names


def _extract_numeric(v: Optional[object]) -> Optional[float]:
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def simplify_record(enriched: Dict[str, object]) -> Dict[str, object]:
    """Flatten a rich merged record to a light-weight view."""
    reaction = enriched.get("reaction", {}) or {}
    substrates = _names_from_compounds(reaction.get("substrates", []) or [])
    products = _names_from_compounds(reaction.get("products", []) or [])

    simplified_measurements: List[Dict[str, object]] = []
    for m in enriched.get("measurements", []) or []:
        kin = m.get("kinetics", {}) or {}
        cond = m.get("conditions", {}) or {}
        simplified_measurements.append(
            {
                "source": m.get("source", {}).get("dataset"),
                "id": m.get("id"),
                "organism": m.get("enzyme", {}).get("organism"),
                "mutant": m.get("enzyme", {}).get("mutant"),
                "direction": m.get("direction"),
                "substrate_ref": m.get("substrate_ref"),
                "substrate_index": m.get("substrate_index"),
                "kcat": _extract_numeric((kin.get("kcat") or {}).get("value")),
                "kcat_unit": (kin.get("kcat") or {}).get("unit"),
                "kcat_std": m.get("kcat_std"),
                "kcat_log10": m.get("kcat_log10"),
                "km": _extract_numeric((kin.get("km") or {}).get("value")),
                "km_unit": (kin.get("km") or {}).get("unit"),
                "km_std_mM": m.get("km_std_mM"),
                "km_log10": m.get("km_log10"),
                "pH": cond.get("pH"),
                "temperature_C": cond.get("temperature"),
                "notes": (m.get("text") or {}).get("descriptor") or (m.get("text") or {}).get(
                    "canonical_description"
                ),
            }
        )

    return {
        "id": enriched.get("id"),
        "primary_ec": enriched.get("primary_ec"),
        "equation": reaction.get("equation_text"),
        "substrates": substrates,
        "products": products,
        "descriptor": (enriched.get("text") or {}).get("descriptor"),
        "measurements": simplified_measurements,
    }


def flatten_record(enriched: Dict[str, object]) -> List[Dict[str, object]]:
    """Produce one row per measurement with reaction context."""
    reaction = enriched.get("reaction", {}) or {}
    substrates = _names_from_compounds(reaction.get("substrates", []) or [])
    products = _names_from_compounds(reaction.get("products", []) or [])
    rows: List[Dict[str, object]] = []
    for m in enriched.get("measurements", []) or []:
        rows.append(
            {
                "reaction_id": enriched.get("id"),
                "primary_ec": enriched.get("primary_ec"),
                "equation": reaction.get("equation_text"),
                "substrates": substrates,
                "products": products,
                "measurement_id": m.get("id"),
                "source": m.get("source", {}).get("dataset"),
                "organism": m.get("enzyme", {}).get("organism"),
                "mutant": m.get("enzyme", {}).get("mutant"),
                "substrate_ref": m.get("substrate_ref"),
                "substrate_index": m.get("substrate_index"),
                "kcat": _extract_numeric((m.get("kinetics") or {}).get("kcat", {}).get("value")),
                "kcat_unit": (m.get("kinetics") or {}).get("kcat", {}).get("unit"),
                "kcat_std": m.get("kcat_std"),
                "kcat_log10": m.get("kcat_log10"),
                "km": _extract_numeric((m.get("kinetics") or {}).get("km", {}).get("value")),
                "km_unit": (m.get("kinetics") or {}).get("km", {}).get("unit"),
                "km_std_mM": m.get("km_std_mM"),
                "km_log10": m.get("km_log10"),
                "pH": (m.get("conditions") or {}).get("pH"),
                "temperature_C": (m.get("conditions") or {}).get("temperature"),
                "direction": m.get("direction"),
                "notes": (m.get("text") or {}).get("descriptor")
                or (m.get("text") or {}).get("canonical_description"),
            }
        )
    return rows


def annotate_measurement(measurement: Dict[str, object], reaction_core: Dict[str, object]) -> Dict[str, object]:
    """Attach substrate_ref, direction, and standardized kinetics to a measurement."""
    m = copy.deepcopy(measurement)
    core_substrates = _names_from_compounds(reaction_core.get("substrates", []) or [])
    meas_substrates = _names_from_compounds((m.get("reaction") or {}).get("substrates", []) or [])
    ref_name, ref_index = _match_substrate(core_substrates, meas_substrates)
    m["substrate_ref"] = ref_name
    m["substrate_index"] = ref_index
    m["direction"] = (m.get("reaction") or {}).get("direction")

    # Standardize kinetics
    kin = m.get("kinetics") or {}
    kcat_val = _extract_numeric((kin.get("kcat") or {}).get("value"))
    kcat_unit = (kin.get("kcat") or {}).get("unit")
    kcat_std = _standardize_kcat(kcat_val, kcat_unit)
    km_val = _extract_numeric((kin.get("km") or {}).get("value"))
    km_unit = (kin.get("km") or {}).get("unit")
    km_std = _standardize_km(km_val, km_unit)

    m["kcat_std"] = kcat_std
    m["kcat_log10"] = math.log10(kcat_std) if kcat_std and kcat_std > 0 else None
    m["km_std_mM"] = km_std
    m["km_log10"] = math.log10(km_std) if km_std and km_std > 0 else None
    return m


def _match_substrate(core_subs: List[str], meas_subs: List[str]) -> tuple[Optional[str], Optional[int]]:
    """Case-insensitive exact matching of measurement substrate to core substrates."""
    core_norm = [s.lower() for s in core_subs]
    for ms in meas_subs:
        if not ms:
            continue
        try_norm = ms.lower()
        if try_norm in core_norm:
            idx = core_norm.index(try_norm)
            return core_subs[idx], idx
    return None, None


def _standardize_kcat(value: Optional[float], unit: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    unit_lower = (unit or "").lower()
    if "min" in unit_lower:
        return value / 60.0
    return value


def _standardize_km(value: Optional[float], unit: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    unit_lower = (unit or "").lower()
    if "um" in unit_lower or "µm" in unit_lower:
        return value / 1000.0
    if unit_lower == "m":
        return value * 1000.0
    # default assume already mM
    return value


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge BKMS reaction cores with EnzyExtract/BRENDA kinetics by EC number into a rich JSONL."
    )
    parser.add_argument("--bkms", type=Path, required=True, help="Path to Reactions_BKMS.csv (TSV).")
    parser.add_argument("--enzyextract", type=Path, help="Path to EnzyExtractDB parquet.")
    parser.add_argument("--brenda", type=Path, help="Path to brenda_kcat_v3 parquet.")
    parser.add_argument("--output", type=Path, required=True, help="Destination JSONL file.")
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on number of BKMS rows to export (for quick sampling).",
    )
    parser.add_argument(
        "--view",
        choices=["full", "minimal", "flat"],
        default="full",
        help="Output format: full (default rich schema), minimal (flattened reaction+measurements), or flat (one row per measurement with reaction context).",
    )
    args = parser.parse_args()

    merged_iter = merge_reactions(
        bkms_path=args.bkms,
        enzy_path=args.enzyextract,
        brenda_path=args.brenda,
        limit=args.limit,
        view=args.view,
    )

    with args.output.open("w", encoding="utf-8") as f:
        for rec in merged_iter:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
