import argparse
import pandas as pd
import re

def load_transcript_ids(transcript_file):
    """Load transcript IDs from a file."""
    with open(transcript_file, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def extract_go_terms(gff_file, transcript_ids):
    """Extract GO terms from GFF file for specified transcript IDs."""
    go_terms_dict = {tid: [] for tid in transcript_ids}
    go_pattern = re.compile(r'em_GOs=([^;]+)')

    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                for tid in transcript_ids:
                    if f"ID={tid}" in line or f"Parent={tid}" in line:
                        match = go_pattern.search(line)
                        if match:
                            go_terms = match.group(1).split(",")
                            go_terms_dict[tid].extend(go_terms)
    return go_terms_dict

def save_to_csv(go_terms_dict, output_file):
    """Save GO terms to CSV."""
    go_df = pd.DataFrame([
        {"Transcript_ID": tid, "GO_Term": go}
        for tid, go_list in go_terms_dict.items()
        for go in set(go_list)
    ])
    go_df.to_csv(output_file, index=False)
    print(f"GO terms extracted and saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract GO terms for specific transcript IDs from a GFF file.")
    parser.add_argument("--transcript_ids", required=True, help="Path to the transcript IDs file (one ID per line).")
    parser.add_argument("--gff", required=True, help="Path to the GFF annotation file.")
    parser.add_argument("--output", required=True, help="Path to save the extracted GO terms CSV file.")

    args = parser.parse_args()

    transcript_ids = load_transcript_ids(args.transcript_ids)
    go_terms_dict = extract_go_terms(args.gff, transcript_ids)
    save_to_csv(go_terms_dict, args.output)

if __name__ == "__main__":
    main()
