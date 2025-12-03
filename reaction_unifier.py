from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional

import numpy as np
import pandas as pd

from unified_reaction_schema import (
    build_base_record,
    build_conditions,
    normalize_list,
    parse_bkms_equation,
    parse_cofactors,
    populate_enzyme_section,
    populate_kinetics,
    populate_reaction,
    populate_text,
)


def _to_native(value: object) -> object:
    if isinstance(value, np.generic):
        return value.item()
    return value


def _clean(value: object) -> Optional[object]:
    """Convert pandas NaN/None/empty strings into None."""
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return None
        if value.size == 1:
            value = value.item()
        else:
            value = value.tolist()
    value = _to_native(value)
    if isinstance(value, float) and pd.isna(value):
        return None
    if isinstance(value, str) and not value.strip():
        return None
    return value


def _substrate_block(
    *,
    name: object,
    name_full: object = None,
    smiles: object = None,
    cid: object = None,
    brenda_id: object = None,
) -> Dict[str, object]:
    return {
        "name": _clean(name),
        "name_full": _clean(name_full),
        "role": "substrate",
        "stoichiometry": 1.0,
        "smiles": _clean(smiles),
        "cid": _clean(cid),
        "brenda_id": _clean(brenda_id),
    }


def _product_block(name: object) -> Dict[str, object]:
    return {
        "name": _clean(name),
        "name_full": None,
        "role": "product",
        "stoichiometry": 1.0,
        "smiles": None,
        "cid": None,
        "brenda_id": None,
    }


def enzyextract_records(path: Path) -> Iterator[Dict[str, object]]:
    df = pd.read_parquet(path)
    for row_index, row in df.iterrows():
        record = build_base_record(
            "EnzyExtractDB",
            int(row_index),
            raw_ids={"brenda_id": _clean(row.get("brenda_id")), "cid": _clean(row.get("cid"))},
        )
        ecs = normalize_list(row.get("enzyme_ecs"))
        record["ec_numbers"] = ecs
        record["primary_ec"] = ecs[0] if ecs else None

        populate_enzyme_section(
            record,
            name=row.get("enzyme"),
            name_full=row.get("enzyme_full"),
            organism=row.get("organism"),
            mutant=row.get("mutant"),
            uniprot_ids=normalize_list(row.get("uniprot")),
            sequence=row.get("sequence"),
            pdb_ids=normalize_list(row.get("pdb")),
        )

        substrates = [
            _substrate_block(
                name=row.get("substrate"),
                name_full=row.get("substrate_full"),
                smiles=row.get("smiles"),
                cid=row.get("cid"),
                brenda_id=row.get("brenda_id"),
            )
        ]
        cofactors = parse_cofactors(row.get("cofactors"))

        populate_reaction(
            record,
            direction="unknown",
            equation_text=None,
            substrates=substrates,
            products=[],
            cofactors=cofactors,
        )

        populate_kinetics(
            record,
            dataset="EnzyExtractDB",
            row_index=int(row_index),
            kcat_value=row.get("kcat_value"),
            kcat_raw=row.get("kcat"),
            km_value=row.get("km_value"),
            km_raw=row.get("km"),
            kcat_km_value=row.get("kcat_km"),
        )

        record["conditions"] = build_conditions(
            dataset="EnzyExtractDB",
            row_index=int(row_index),
            pH=row.get("pH"),
            temperature=row.get("temperature"),
            temperature_unit="C",
            comments=row.get("solution"),
        )

        populate_text(
            record,
            descriptor=row.get("descriptor"),
            canonical_description=row.get("canonical"),
            notes=row.get("other"),
        )

        yield record


def brenda_records(path: Path) -> Iterator[Dict[str, object]]:
    df = pd.read_parquet(path)
    for row_index, row in df.iterrows():
        record = build_base_record(
            "BRENDA",
            int(row_index),
            raw_ids={"pmid": _clean(row.get("pmid")), "protein_id": _clean(row.get("protein_id"))},
        )
        ecs = normalize_list(row.get("ec"))
        record["ec_numbers"] = ecs or [row.get("ec")] if _clean(row.get("ec")) else []
        record["primary_ec"] = record["ec_numbers"][0] if record["ec_numbers"] else None

        populate_enzyme_section(
            record,
            name=row.get("enzyme"),
            organism=row.get("organism_name"),
            mutant=row.get("mutant"),
            uniprot_ids=normalize_list(row.get("accessions")),
        )

        substrates = [
            _substrate_block(
                name=row.get("substrate"),
            )
        ]

        populate_reaction(
            record,
            direction="unknown",
            equation_text=None,
            substrates=substrates,
            products=[],
            cofactors=[],
        )

        populate_kinetics(
            record,
            dataset="BRENDA",
            row_index=int(row_index),
            kcat_value=row.get("turnover_number"),
            kcat_raw=row.get("turnover_number"),
            km_value=row.get("km_value"),
            km_raw=row.get("km_value"),
            kcat_km_value=row.get("kcat_km"),
        )

        record["conditions"] = build_conditions(
            dataset="BRENDA",
            row_index=int(row_index),
            pH=row.get("pH"),
            temperature=row.get("temperature"),
            temperature_unit="C",
            comments=row.get("comments"),
        )

        populate_text(
            record,
            descriptor=row.get("comments"),
            canonical_description=None,
            notes=row.get("organism_comments"),
        )

        yield record


def bkms_records(path: Path) -> Iterator[Dict[str, object]]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    for row_index, row in df.iterrows():
        record = build_base_record(
            "BKMS",
            int(row_index),
            raw_ids={
                "bkms_reaction_id": _clean(row.get("Reaction_ID_BRENDA")),
                "kegg_id": _clean(row.get("Reaction_ID_KEGG")),
                "metacyc_id": _clean(row.get("Reaction_ID_MetaCyc")),
                "sabiork_id": _clean(row.get("Reaction_ID_SABIO_RK")),
            },
        )
        ecs = normalize_list(row.get("EC_Number"))
        record["ec_numbers"] = ecs
        record["primary_ec"] = ecs[0] if ecs else None

        populate_enzyme_section(record, name=row.get("Recommended_Name"))

        reaction_parts = parse_bkms_equation(row.get("Reaction"))
        populate_reaction(
            record,
            direction=reaction_parts["direction"],
            equation_text=reaction_parts["equation_text"],
            substrates=reaction_parts["substrates"],
            products=reaction_parts["products"],
            cofactors=[],
        )

        record["conditions"] = build_conditions(
            dataset="BKMS",
            row_index=int(row_index),
            temperature=None,
            pH=None,
            temperature_unit=None,
            comments=None,
        )

        populate_text(
            record,
            descriptor=row.get("Remark"),
            canonical_description=None,
            notes=row.get("Commentary_KEGG") or row.get("Commentary_MetaCyc"),
        )

        yield record


def _writer(records: Iterable[Dict[str, object]], output: Path) -> None:
    with output.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Project EnzyExtract/BRENDA/BKMS datasets into a unified LLM-friendly JSONL schema."
    )
    parser.add_argument(
        "--dataset",
        choices=["enzyextract", "brenda", "bkms"],
        required=True,
        help="Which source dataset to convert.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input file path (parquet for EnzyExtract/BRENDA, TSV for BKMS).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination JSONL file path.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit to first N rows for quick sampling.",
    )
    args = parser.parse_args()

    if args.dataset == "enzyextract":
        generator: Iterable[Dict[str, object]] = enzyextract_records(args.input)
    elif args.dataset == "brenda":
        generator = brenda_records(args.input)
    else:
        generator = bkms_records(args.input)

    if args.limit is not None:
        def limited(iterable: Iterable[Dict[str, object]], n: int) -> Iterator[Dict[str, object]]:
            count = 0
            for item in iterable:
                if count >= n:
                    break
                count += 1
                yield item

        generator = limited(generator, args.limit)

    _writer(generator, args.output)


if __name__ == "__main__":
    main()
