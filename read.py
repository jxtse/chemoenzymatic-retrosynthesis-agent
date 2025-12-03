import gzip
import json
import pandas as pd
import numpy as np
from pathlib import Path
import os
import sys

# Set encoding for Windows environment
if sys.platform == 'win32':
    os.environ['PYTHONIOENCODING'] = 'utf-8'


def safe_print(text):
    """Print text safely, handling encoding issues"""
    try:
        print(text)
    except (UnicodeEncodeError, UnicodeDecodeError):
        # Remove problematic characters
        safe_text = str(text).encode('utf-8', errors='ignore').decode('utf-8', errors='ignore')
        print(safe_text)


def read_bkms_json(file_path="bkms-data/reactions/bkms-reactions.json.gz", num_examples=5):
    """Read BKMS JSON data from gzipped file"""
    print(f"\n=== BKMS JSON Data ===")
    try:
        with gzip.open(file_path, "rt", encoding="utf8") as f:
            data = json.load(f)

        print(f"Total entries: {len(data)}")
        print(f"Showing first {num_examples} examples:")

        # Check if data is a list of dictionaries or a single dictionary
        if isinstance(data, dict) and 'data' in data:
            data = data['data']  # Extract actual data if nested

        for i, ex in enumerate(data[:num_examples]):
            if isinstance(ex, dict):
                print(f"\nExample {i+1}:")
                print(f"  ID: {ex.get('id', ex.get('ID', 'N/A'))}")
                print(f"  Name: {ex.get('name', ex.get('Name', 'N/A'))}")
                print(f"  EC: {ex.get('ec', ex.get('EC', 'N/A'))}")
                print(f"  Reaction: {ex.get('reaction', ex.get('Reaction', 'N/A'))}")
            else:
                print(f"\nExample {i+1}: {ex}")

        return data
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading BKMS JSON: {e}")
        return None


def read_bkms_csv(file_path="Reactions_BKMS.csv", num_examples=5):
    """Read BKMS CSV data"""
    print(f"\n=== BKMS CSV Data ===")
    try:
        # Fix the mixed types warning by specifying low_memory=False
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        print(f"Total reactions: {len(df)}")
        print(f"Columns: {list(df.columns)}")
        print(f"Showing first {num_examples} examples:")

        for i, (_, row) in enumerate(df.head(num_examples).iterrows()):
            print(f"\nExample {i+1}:")
            print(f"  ID: {row.get('ID', 'N/A')}")
            print(f"  EC Number: {row.get('EC_Number', 'N/A')}")
            print(f"  Name: {row.get('Recommended_Name', 'N/A')}")
            print(f"  Reaction: {row.get('Reaction', 'N/A')}")
            print(f"  BRENDA ID: {row.get('Reaction_ID_BRENDA', 'N/A')}")
            print(f"  KEGG ID: {row.get('Reaction_ID_KEGG', 'N/A')}")

        return df
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading BKMS CSV: {e}")
        return None


def read_brenda_parquet(file_path="brenda_kcat_v3.parquet", num_examples=5):
    """Read BRENDA database parquet data"""
    print(f"\n=== BRENDA Database (Parquet) ===")
    try:
        df = pd.read_parquet(file_path)
        print(f"Total entries: {len(df)}")
        print(f"Columns: {list(df.columns)}")
        print(f"Data types:\n{df.dtypes}")
        print(f"Showing first {num_examples} examples:")

        for i, (_, row) in enumerate(df.head(num_examples).iterrows()):
            print(f"\nExample {i+1}:")
            for col in df.columns:
                value = row[col]
                if pd.isna(value):
                    value = "N/A"
                elif isinstance(value, str):
                    try:
                        if len(value) > 100:
                            value = value[:100] + "..."
                        # Test if the string can be encoded/decoded safely
                        value.encode('utf-8').decode('utf-8')
                    except (UnicodeEncodeError, UnicodeDecodeError):
                        # Replace problematic characters
                        value = value.encode('utf-8', errors='ignore').decode('utf-8', errors='ignore')
                        if len(value) > 100:
                            value = value[:100] + "..."
                print(f"  {col}: {value}")

        return df
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading BRENDA parquet: {e}")
        return None


def read_enzyextract_parquet(file_path="EnzyExtractDB_176463.parquet", num_examples=5):
    """Read EnzyExtractDB parquet data"""
    safe_print(f"\n=== EnzyExtractDB Database (Parquet) ===")
    try:
        df = pd.read_parquet(file_path)
        safe_print(f"Total entries: {len(df)}")
        safe_print(f"Columns: {list(df.columns)}")
        safe_print(f"Data types:\n{df.dtypes}")
        safe_print(f"Showing first {num_examples} examples:")

        examples_data = []
        for i, (_, row) in enumerate(df.head(num_examples).iterrows()):
            try:
                safe_print(f"\nExample {i+1}:")
                example_data = {}
                for col in df.columns:
                    try:
                        value = row[col]
                        if pd.isna(value):
                            value = "N/A"
                        elif isinstance(value, str):
                            # Handle encoding issues by replacing problematic characters
                            try:
                                if len(value) > 100:
                                    value = value[:100] + "..."
                                # Test if the string can be encoded/decoded safely
                                value.encode('utf-8').decode('utf-8')
                            except (UnicodeEncodeError, UnicodeDecodeError):
                                # Replace problematic characters
                                value = value.encode('utf-8', errors='ignore').decode('utf-8', errors='ignore')
                                if len(value) > 100:
                                    value = value[:100] + "..."
                        safe_print(f"  {col}: {value}")
                        example_data[col] = value
                    except Exception as col_error:
                        safe_print(f"  {col}: [Encoding error - {col_error}]")
                        example_data[col] = f"[Encoding error - {col_error}]"
                examples_data.append(example_data)
            except Exception as row_error:
                safe_print(f"Example {i+1}: [Row processing error - {row_error}]")
                continue

        return df, examples_data
    except FileNotFoundError:
        safe_print(f"File not found: {file_path}")
        return None, None
    except Exception as e:
        safe_print(f"Error reading EnzyExtractDB parquet: {e}")
        return None, None


def save_bkms_json_examples(data, file_path="bkms_json_examples.txt", num_examples=10):
    """Save BKMS JSON examples to local file"""
    safe_print(f"\n=== Saving BKMS JSON Examples ===")
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f"BKMS JSON Examples (First {min(num_examples, len(data))} entries)\n")
            f.write("=" * 50 + "\n\n")

            for i, ex in enumerate(data[:num_examples]):
                f.write(f"Example {i+1}:\n")
                if isinstance(ex, dict):
                    f.write(f"  ID: {ex.get('id', ex.get('ID', 'N/A'))}\n")
                    f.write(f"  Name: {ex.get('name', ex.get('Name', 'N/A'))}\n")
                    f.write(f"  EC: {ex.get('ec', ex.get('EC', 'N/A'))}\n")
                    f.write(f"  Reaction: {ex.get('reaction', ex.get('Reaction', 'N/A'))}\n")
                else:
                    f.write(f"  Data: {ex}\n")
                f.write("\n")

        safe_print(f"BKMS JSON examples saved to: {file_path}")
        return True
    except Exception as e:
        safe_print(f"Error saving BKMS JSON examples: {e}")
        return False


def save_bkms_csv_examples(df, file_path="bkms_csv_examples.txt", num_examples=10):
    """Save BKMS CSV examples to local file"""
    safe_print(f"\n=== Saving BKMS CSV Examples ===")
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f"BKMS CSV Examples (First {min(num_examples, len(df))} entries)\n")
            f.write("=" * 50 + "\n\n")

            for i, (_, row) in enumerate(df.head(num_examples).iterrows()):
                f.write(f"Example {i+1}:\n")
                f.write(f"  ID: {row.get('ID', 'N/A')}\n")
                f.write(f"  EC Number: {row.get('EC_Number', 'N/A')}\n")
                f.write(f"  Name: {row.get('Recommended_Name', 'N/A')}\n")
                f.write(f"  Reaction: {row.get('Reaction', 'N/A')}\n")
                f.write(f"  BRENDA ID: {row.get('Reaction_ID_BRENDA', 'N/A')}\n")
                f.write(f"  KEGG ID: {row.get('Reaction_ID_KEGG', 'N/A')}\n")
                f.write(f"  KEGG Pathway: {row.get('KEGG_Pathway_Name', 'N/A')}\n")
                f.write(f"  MetaCyc Pathway: {row.get('MetaCyc_Pathway_Name', 'N/A')}\n")
                f.write("\n")

        safe_print(f"BKMS CSV examples saved to: {file_path}")
        return True
    except Exception as e:
        safe_print(f"Error saving BKMS CSV examples: {e}")
        return False


def save_brenda_examples(df, file_path="brenda_parquet_examples.txt", num_examples=10):
    """Save BRENDA parquet examples to local file"""
    safe_print(f"\n=== Saving BRENDA Parquet Examples ===")
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f"BRENDA Parquet Examples (First {min(num_examples, len(df))} entries)\n")
            f.write("=" * 50 + "\n\n")

            for i, (_, row) in enumerate(df.head(num_examples).iterrows()):
                f.write(f"Example {i+1}:\n")
                for col in df.columns:
                    value = row[col]
                    if pd.isna(value):
                        value = "N/A"
                    elif isinstance(value, str) and len(value) > 100:
                        value = value[:100] + "..."
                    f.write(f"  {col}: {value}\n")
                f.write("\n")

        safe_print(f"BRENDA parquet examples saved to: {file_path}")
        return True
    except Exception as e:
        safe_print(f"Error saving BRENDA parquet examples: {e}")
        return False


def save_enzyextract_examples(examples_data, file_path="enzyextractdb_examples.txt", num_examples=10):
    """Save EnzyExtractDB parquet examples to local file"""
    safe_print(f"\n=== Saving EnzyExtractDB Examples ===")
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f"EnzyExtractDB Parquet Examples (First {min(num_examples, len(examples_data))} entries)\n")
            f.write("=" * 50 + "\n\n")

            for i, example in enumerate(examples_data[:num_examples]):
                f.write(f"Example {i+1}:\n")
                for col, value in example.items():
                    f.write(f"  {col}: {value}\n")
                f.write("\n")

        safe_print(f"EnzyExtractDB parquet examples saved to: {file_path}")
        return True
    except Exception as e:
        safe_print(f"Error saving EnzyExtractDB parquet examples: {e}")
        return False


def compare_datasets(bkms_json_data, bkms_csv_data, brenda_data):
    """Compare different datasets"""
    print(f"\n=== Dataset Comparison ===")

    # Basic statistics
    print(f"\nDataset Sizes:")
    if bkms_json_data is not None:
        print(f"  BKMS JSON: {len(bkms_json_data)} entries")
    if bkms_csv_data is not None:
        print(f"  BKMS CSV: {len(bkms_csv_data)} reactions")
    if brenda_data is not None:
        print(f"  BRENDA: {len(brenda_data)} entries")

    # EC number comparison
    if bkms_csv_data is not None and brenda_data is not None:
        print(f"\nEC Number Analysis:")

        # BKMS EC numbers
        bkms_ecs = set()
        for ec in bkms_csv_data['EC_Number'].dropna():
            if pd.notna(ec):
                bkms_ecs.add(str(ec).strip())

        # BRENDA EC numbers (assuming there's an EC column)
        brenda_ec_cols = [col for col in brenda_data.columns if 'ec' in col.lower()]
        if brenda_ec_cols:
            brenda_ecs = set()
            for col in brenda_ec_cols:
                for ec in brenda_data[col].dropna():
                    if pd.notna(ec):
                        brenda_ecs.add(str(ec).strip())

            print(f"  BKMS unique EC numbers: {len(bkms_ecs)}")
            print(f"  BRENDA unique EC numbers: {len(brenda_ecs)}")
            print(f"  Common EC numbers: {len(bkms_ecs.intersection(brenda_ecs))}")
            print(f"  BKMS-only EC numbers: {len(bkms_ecs - brenda_ecs)}")
            print(f"  BRENDA-only EC numbers: {len(brenda_ecs - bkms_ecs)}")

            # Show some examples of common and unique EC numbers
            common_ecs = list(bkms_ecs.intersection(brenda_ecs))[:5]
            bkms_only_ecs = list(bkms_ecs - brenda_ecs)[:5]
            brenda_only_ecs = list(brenda_ecs - bkms_ecs)[:5]

            print(f"\n  Common EC numbers (examples): {common_ecs}")
            print(f"  BKMS-only EC numbers (examples): {bkms_only_ecs}")
            print(f"  BRENDA-only EC numbers (examples): {brenda_only_ecs}")

    # Reaction type analysis
    if bkms_csv_data is not None:
        print(f"\nBKMS Reaction Types:")
        ec_types = bkms_csv_data['EC_Number'].dropna()
        if len(ec_types) > 0:
            main_classes = [str(ec).split('.')[0] for ec in ec_types if str(ec).count('.') >= 1]
            if main_classes:
                from collections import Counter
                class_counts = Counter(main_classes)
                print(f"  Main EC classes: {dict(class_counts)}")

    print(f"\nData Quality Notes:")
    if bkms_csv_data is not None:
        missing_reactions = bkms_csv_data['Reaction'].isna().sum()
        print(f"  BKMS CSV missing reactions: {missing_reactions}/{len(bkms_csv_data)}")

    if brenda_data is not None:
        for col in brenda_data.columns:
            missing = brenda_data[col].isna().sum()
            print(f"  BRENDA {col} missing: {missing}/{len(brenda_data)}")


if __name__ == "__main__":
    # Set number of examples to 10 for saving
    NUM_EXAMPLES = 10

    # Read all datasets with 10 examples
    bkms_json_data = read_bkms_json(num_examples=NUM_EXAMPLES)
    bkms_csv_data = read_bkms_csv(num_examples=NUM_EXAMPLES)
    brenda_data = read_brenda_parquet(num_examples=NUM_EXAMPLES)
    enzyextract_df, enzyextract_data = read_enzyextract_parquet(num_examples=NUM_EXAMPLES)

    # Save examples to local files
    if bkms_json_data is not None:
        save_bkms_json_examples(bkms_json_data, num_examples=NUM_EXAMPLES)

    if bkms_csv_data is not None:
        save_bkms_csv_examples(bkms_csv_data, num_examples=NUM_EXAMPLES)

    if brenda_data is not None:
        save_brenda_examples(brenda_data, num_examples=NUM_EXAMPLES)

    if enzyextract_data is not None:
        save_enzyextract_examples(enzyextract_data, num_examples=NUM_EXAMPLES)

    # Compare datasets (only compare the original three)
    compare_datasets(bkms_json_data, bkms_csv_data, brenda_data)

    safe_print(f"\n=== Summary ===")
    safe_print(f"Saved {NUM_EXAMPLES} examples from each dataset to local text files:")
    safe_print(f"  - bkms_json_examples.txt")
    safe_print(f"  - bkms_csv_examples.txt")
    safe_print(f"  - brenda_parquet_examples.txt")
    safe_print(f"  - enzyextractdb_examples.txt")
