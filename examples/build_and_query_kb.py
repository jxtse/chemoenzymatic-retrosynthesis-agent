#!/usr/bin/env python3
"""
Example: Build and Query the Unified Knowledge Base

This script demonstrates:
1. Building the knowledge base from multiple sources
2. Querying by EC number, compound, and reaction
3. Extracting kinetic parameters
4. Getting statistics
"""

from pathlib import Path

from knowledge_base import KnowledgeBaseAPI, KnowledgeBaseBuilder, KnowledgeBaseConfig


def build_knowledge_base():
    """Build the knowledge base from configured sources."""
    print("=" * 70)
    print("STEP 1: Building Knowledge Base")
    print("=" * 70)

    # Create configuration
    config = KnowledgeBaseConfig.default()

    # Customize paths (adjust to your local files)
    config.bkms.path = Path("Reactions_BKMS.csv")
    config.brenda.path = Path("brenda_kcat_v3.parquet")
    config.enzyextract.path = Path("EnzyExtractDB_176463.parquet")
    config.retrobiocat.path = Path("RetroBioCat/rxns_yaml.yaml")

    # Set output
    config.output_dir = Path("kb_output")
    config.output_format = "jsonl"
    config.compress_output = False

    # Build
    builder = KnowledgeBaseBuilder(config)

    # Build without API enrichment (faster)
    output_file = builder.build_full_knowledge_base(enrich=False)

    print(f"\n✓ Knowledge base built: {output_file}\n")
    return output_file


def query_examples(kb_path: Path):
    """Demonstrate various query types."""
    print("=" * 70)
    print("STEP 2: Querying Knowledge Base")
    print("=" * 70)

    # Load knowledge base
    api = KnowledgeBaseAPI(kb_path)
    api.load()

    # Example 1: Query by EC number
    print("\n[Example 1] Query by EC number: 1.1.1.1")
    print("-" * 70)
    result = api.query_by_ec("1.1.1.1")
    print(f"Found {result['count']} records for EC 1.1.1.1")

    if result['records']:
        record = result['records'][0]
        print(f"\nFirst record:")
        print(f"  ID: {record['id']}")
        print(f"  Source: {record['source']['dataset']}")
        print(f"  Enzyme: {record.get('enzyme', {}).get('name', 'N/A')}")
        print(f"  Reaction: {record.get('reaction', {}).get('equation_text', 'N/A')}")

    # Example 2: Query by compound
    print("\n\n[Example 2] Query by compound: glucose")
    print("-" * 70)
    result = api.query_by_compound("glucose", role="substrate")
    print(f"Found {result['count']} reactions with glucose as substrate")

    if result['records']:
        for i, record in enumerate(result['records'][:3]):  # Show first 3
            print(f"\n  {i+1}. {record.get('reaction', {}).get('equation_text', 'N/A')}")
            print(f"     EC: {record.get('primary_ec', 'N/A')}")

    # Example 3: Get kinetics
    print("\n\n[Example 3] Get kinetic parameters for EC 1.1.1.1")
    print("-" * 70)
    kinetics = api.get_kinetics(ec_number="1.1.1.1")
    print(f"Found {kinetics['count']} kinetic measurements")

    if kinetics['kinetics']:
        for i, k in enumerate(kinetics['kinetics'][:5]):  # Show first 5
            print(f"\n  {i+1}. {k['substrate']} ({k.get('organism', 'N/A')})")
            if k['kcat']:
                print(f"     kcat: {k['kcat']} {k['kcat_unit']}")
            if k['km']:
                print(f"     Km:   {k['km']} {k['km_unit']}")
            print(f"     Source: {k['source']}")

    # Example 4: Query by reaction
    print("\n\n[Example 4] Query by reaction: ATP + glucose -> ADP + glucose-6-phosphate")
    print("-" * 70)
    result = api.query_by_reaction(
        substrates=["ATP", "glucose"],
        products=["ADP", "glucose-6-phosphate"],
    )
    print(f"Found {result['count']} matching reactions")

    if result['records']:
        for i, record in enumerate(result['records'][:3]):
            print(f"\n  {i+1}. EC {record.get('primary_ec', 'N/A')}")
            print(f"     {record.get('reaction', {}).get('equation_text', 'N/A')}")

    # Example 5: Statistics
    print("\n\n[Example 5] Knowledge Base Statistics")
    print("-" * 70)
    stats = api.get_statistics()
    print(f"Total records:          {stats['total_records']}")
    print(f"Unique EC numbers:      {stats['unique_ec_numbers']}")
    print(f"Unique compounds:       {stats['unique_compounds']}")
    print(f"Unique organisms:       {stats['unique_organisms']}")
    print(f"Records with kinetics:  {stats['records_with_kinetics']}")
    print(f"Records with sequence:  {stats['records_with_sequence']}")

    print("\nRecords by source:")
    for source, count in stats['by_source'].items():
        print(f"  {source:20s}: {count:6d}")

    print("\nRecords by EC class:")
    for ec_class, count in stats['by_ec_class'].items():
        print(f"  {ec_class:20s}: {count:6d}")

    # Export example query result
    print("\n\n[Example 6] Export query result to JSON")
    print("-" * 70)
    output_path = Path("kb_output/query_result_ec_1.1.1.1.json")
    result = api.query_by_ec("1.1.1.1")
    api.export_to_json(output_path, result)
    print(f"✓ Exported to {output_path}")


def main():
    """Main execution."""
    print("\n" + "=" * 70)
    print("Unified Chemical-Enzyme-Kinetics Knowledge Base Example")
    print("=" * 70)

    # Build knowledge base
    kb_path = build_knowledge_base()

    # Query examples
    query_examples(kb_path)

    print("\n" + "=" * 70)
    print("Example completed!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
