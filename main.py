import argparse

from brenda_client import get_km_values, get_reactions, load_env_from_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch Km values or reactions from the BRENDA SOAP API."
    )
    parser.add_argument(
        "--mode",
        choices=["km", "reaction"],
        default="km",
        help="Choose 'km' (default) to fetch Km values or 'reaction' to fetch reaction equations.",
    )
    parser.add_argument(
        "--ec",
        required=True,
        help="EC number to query, for example 1.1.1.1",
    )
    parser.add_argument(
        "--organism",
        default="*",
        help="Organism filter; use '*' to match all organisms (default)",
    )
    parser.add_argument(
        "--substrate",
        default="*",
        help="Substrate filter; use '*' to match all substrates (default)",
    )
    parser.add_argument(
        "--reaction",
        default="*",
        help="Reaction filter for reaction mode; use '*' to match all reactions (default)",
    )
    args = parser.parse_args()

    load_env_from_file()
    if args.mode == "km":
        entries = get_km_values(
            ec_number=args.ec, organism=args.organism, substrate=args.substrate
        )
    else:
        entries = get_reactions(
            ec_number=args.ec, organism=args.organism, reaction=args.reaction
        )

    if not entries:
        print("No entries returned for the given filters.")
        return

    for entry in entries:
        print(entry)


if __name__ == "__main__":
    main()
