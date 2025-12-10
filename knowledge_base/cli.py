"""Command-line interface for knowledge base operations."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .api import KnowledgeBaseAPI
from .builder import KnowledgeBaseBuilder
from .config import KnowledgeBaseConfig

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def cmd_build(args: argparse.Namespace) -> int:
    """Build knowledge base from data sources."""
    # Load config
    if args.config:
        config = KnowledgeBaseConfig.from_yaml(Path(args.config))
    else:
        config = KnowledgeBaseConfig.default()

    # Override output settings
    if args.output:
        config.output_dir = Path(args.output)
    if args.format:
        config.output_format = args.format

    # Build
    try:
        builder = KnowledgeBaseBuilder(config)
        output_file = builder.build_full_knowledge_base(enrich=args.enrich)
        print(f"\n✓ Knowledge base built successfully: {output_file}")
        return 0
    except Exception as e:
        logger.error(f"Build failed: {e}", exc_info=True)
        return 1


def cmd_query(args: argparse.Namespace) -> int:
    """Query the knowledge base."""
    kb_path = Path(args.kb)
    if not kb_path.exists():
        print(f"Error: Knowledge base not found: {kb_path}")
        return 1

    try:
        api = KnowledgeBaseAPI(kb_path)
        api.load()

        # Execute query
        if args.ec:
            result = api.query_by_ec(args.ec)
        elif args.compound:
            result = api.query_by_compound(args.compound, role=args.role)
        elif args.reaction:
            substrates = args.substrates.split(",") if args.substrates else None
            products = args.products.split(",") if args.products else None
            result = api.query_by_reaction(substrates, products)
        elif args.kinetics:
            result = api.get_kinetics(
                ec_number=args.ec_filter,
                substrate=args.substrate_filter,
                organism=args.organism_filter,
            )
        elif args.stats:
            result = api.get_statistics()
        else:
            print("Error: No query specified. Use --ec, --compound, --reaction, --kinetics, or --stats")
            return 1

        # Output result
        if args.output:
            api.export_to_json(Path(args.output), result, indent=2)
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))

        return 0

    except Exception as e:
        logger.error(f"Query failed: {e}", exc_info=True)
        return 1


def cmd_init_config(args: argparse.Namespace) -> int:
    """Initialize default configuration file."""
    output_path = Path(args.output)

    if output_path.exists() and not args.force:
        print(f"Error: Config file already exists: {output_path}")
        print("Use --force to overwrite")
        return 1

    try:
        config = KnowledgeBaseConfig.default()
        config.save_yaml(output_path)
        print(f"✓ Created default config: {output_path}")
        print("\nEdit this file to customize your data sources and settings.")
        return 0
    except Exception as e:
        logger.error(f"Failed to create config: {e}")
        return 1


def main() -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Unified Chemical-Enzyme-Kinetics Knowledge Base",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # Build command
    build_parser = subparsers.add_parser(
        "build",
        help="Build knowledge base from data sources",
    )
    build_parser.add_argument(
        "--config",
        type=str,
        help="Configuration YAML file (uses defaults if not specified)",
    )
    build_parser.add_argument(
        "--output",
        type=str,
        help="Output directory (overrides config)",
    )
    build_parser.add_argument(
        "--format",
        choices=["jsonl", "parquet"],
        help="Output format (overrides config)",
    )
    build_parser.add_argument(
        "--enrich",
        action="store_true",
        help="Enrich with API data (KEGG, UniProt, PubChem)",
    )
    build_parser.set_defaults(func=cmd_build)

    # Query command
    query_parser = subparsers.add_parser(
        "query",
        help="Query the knowledge base",
    )
    query_parser.add_argument(
        "--kb",
        required=True,
        help="Path to knowledge base file (JSONL or Parquet)",
    )

    # Query types
    query_group = query_parser.add_mutually_exclusive_group()
    query_group.add_argument("--ec", help="Query by EC number")
    query_group.add_argument("--compound", help="Query by compound name")
    query_group.add_argument("--reaction", action="store_true", help="Query by reaction")
    query_group.add_argument("--kinetics", action="store_true", help="Get kinetic parameters")
    query_group.add_argument("--stats", action="store_true", help="Get knowledge base statistics")

    # Query filters
    query_parser.add_argument("--role", choices=["substrate", "product"], help="Filter by compound role")
    query_parser.add_argument("--substrates", help="Comma-separated substrate names")
    query_parser.add_argument("--products", help="Comma-separated product names")
    query_parser.add_argument("--ec-filter", help="Filter kinetics by EC number")
    query_parser.add_argument("--substrate-filter", help="Filter kinetics by substrate")
    query_parser.add_argument("--organism-filter", help="Filter kinetics by organism")

    # Output
    query_parser.add_argument(
        "--output",
        help="Output JSON file (prints to stdout if not specified)",
    )
    query_parser.set_defaults(func=cmd_query)

    # Init config command
    init_parser = subparsers.add_parser(
        "init-config",
        help="Initialize default configuration file",
    )
    init_parser.add_argument(
        "--output",
        default="kb_config.yaml",
        help="Output config file path (default: kb_config.yaml)",
    )
    init_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing config file",
    )
    init_parser.set_defaults(func=cmd_init_config)

    # Parse and execute
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
