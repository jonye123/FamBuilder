import os
import argparse
import pandas as pd
import shutil
import statistics
from Bio import SeqIO

# Check that orthofinder_to_fasta is available
try:
    import orthofinder_to_fasta
except ImportError:
    print("Warning: orthofinder_to_fasta.py not found. Iteration 1 extraction might fail if not handled elsewhere.")

def get_args():
    parser = argparse.ArgumentParser(description="Filter HOGs and generate stats.")
    parser.add_argument("--input", required=True, help="Path to N0.tsv (or input file)")
    parser.add_argument("--fasta_dir", help="Directory containing original proteome FASTA files")
    parser.add_argument("--project_name", required=True, help="Name of the project (e.g., FungiFam, GrassFam)")
    parser.add_argument("--iteration", type=int, required=True, help="Current iteration number")
    parser.add_argument("--min_species", type=int, default=4, help="Minimum number of species required to keep a HOG")
    parser.add_argument("--grass_filter", action="store_true", help="Apply strict BOP/PACMAD coverage rules")
    return parser.parse_args()

def calculate_stats(hog_directory, output_tsv):
    """
    Scans a directory of FASTA files and generates a statistics table.
    Calculates: Min/Max/Mean/Median Length, Species Count, Gene Count.
    """
    if not os.path.exists(hog_directory):
        print(f"Directory {hog_directory} does not exist. Creating empty stats.")
        return pd.DataFrame()
    
    # Detect species from the first few files to build headers
    all_species = set()
    sample_files = [f for f in os.listdir(hog_directory) if f.endswith(".fasta")][:100] # take first 100 fasta files

    for filename in sample_files:
        try:
            for record in SeqIO.parse(os.path.join(hog_directory, filename), "fasta"):
                all_species.add(record.id.split('|')[0])
        except:
            pass
    
    species_list = sorted(list(all_species))
    
    # Prepare output file
    with open(output_tsv, 'w') as f_out:
        species_columns = "\t".join(species_list)
        f_out.write(f"HOG\tmin\tmax\tmedian\tmean\tsd\tspecies_count\tgene_count\tlistOfLengths\t{species_columns}\n")

        # Iterate through every FASTA file
        for filename in os.listdir(hog_directory):
            if filename.endswith(".fasta"):
                # Extract HOG ID (Handling generic N0.HOG001.fasta or HOG001.fasta)
                hog_id = filename.replace(".fasta", "").replace("N0.", "")
                hog_file = os.path.join(hog_directory, filename)

                species = set()
                lengths = []
                
                try:
                    for record in SeqIO.parse(hog_file, "fasta"):
                        lengths.append(len(record.seq))
                        # robust split in case ID format varies
                        parts = record.id.split('|')
                        if len(parts) > 0:
                            species.add(parts[0])
                except Exception as e:
                    print(f"Skipping corrupt file {filename}: {e}")
                    continue

                if not lengths: 
                    continue

                # Calculate Stats
                gene_count = len(lengths)
                min_len = min(lengths)
                max_len = max(lengths)
                med_len = statistics.median(lengths)
                mean_len = statistics.mean(lengths)
                sd_len = statistics.stdev(lengths) if gene_count > 1 else 0
                
                lengths_str = ','.join(map(str, lengths))
                species_count = len(species)
                
                # Binary presence/absence for species columns
                species_binary = [str(1 if s in species else 0) for s in species_list]

                # Write row
                f_out.write(f"{hog_id}\t{min_len}\t{max_len}\t{med_len}\t{mean_len}\t{sd_len}\t{species_count}\t{gene_count}\t{lengths_str}\t{'\t'.join(species_binary)}\n")
    
    return pd.read_csv(output_tsv, sep='\t')

def clean_genes_from_hogs(hog_directory, genes_to_remove):
    """
    Opens every FASTA file in hog_directory.
    Removes sequences whose IDs are in the 'genes_to_remove' set.
    """
    print(f"Cleaning {len(genes_to_remove)} bad genes from FASTA files...")
    
    count_removed = 0
    for filename in os.listdir(hog_directory):
        if not filename.endswith(".fasta"): continue
        
        file_path = os.path.join(hog_directory, filename)
        temp_path = file_path + ".tmp"
        
        modified = False
        kept_sequences = []
        
        # Read and Filter
        for record in SeqIO.parse(file_path, "fasta"):
            # Robust ID extraction: Handle 'Species|Gene' vs 'Gene'
            gene_id = record.id.split('|')[-1] 
            full_id = record.id 

            # Check if this gene is in our removal list
            if gene_id in genes_to_remove or full_id in genes_to_remove:
                modified = True
                count_removed += 1
            else:
                kept_sequences.append(record)
        
        # Write back only if we changed something
        if modified:
            if len(kept_sequences) > 0:
                SeqIO.write(kept_sequences, temp_path, "fasta")
                os.replace(temp_path, file_path)
            else:
                # If we deleted ALL genes, delete the file
                os.remove(file_path)

    print(f"Removed {count_removed} sequences.")

def main():
    args = get_args()
    
    # Define grass clades from original GrassFam
    pacmad = ['Svi', 'Pvi', 'Pha', 'Msi', 'Sbi', 'Zma', 'Ecu', 'Pau']
    bop = ['Osa', 'Ogl', 'Ola', 'Lpe', 'Hvu', 'Ata', 'Tae', 'Bdi']

    print(f"Filtering HOGs (GrassFam Mode: {args.grass_filter})...")
    
    # --- ITERATION 1: INITIAL SETUP ---
    if args.iteration == 1:
        print(f"Iteration 1: Setting up {args.project_name}...")
        
        # In the notebook, this step usually extracts sequences from N0.tsv
        # We call the helper script 'orthofinder_to_fasta'
        if os.path.exists(args.input):
            print(f"Parsing {args.input}...")
            
            # Use the imported module to parse and extract
            # Note: You might need to adjust this function call depending on 
            # exactly how your orthofinder_to_fasta.py is written.
            try:
                orthogroups = orthofinder_to_fasta.parse_orthogroups(args.input)
                # Output to generic 'hogDirectory'
                orthofinder_to_fasta.extract_sequences(orthogroups, args.fasta_dir, "./hogDirectory")
            except AttributeError:
                print("Error: Could not call methods from orthofinder_to_fasta. Make sure the file exists and has 'parse_orthogroups' and 'extract_sequences' functions.")
                # Fallback: Create directory if it doesn't exist so stats don't crash
                os.makedirs("hogDirectory", exist_ok=True)
        else:
            print(f"Error: Input file {args.input} not found.")

    # --- ITERATION > 1: REFINEMENT ---
    else:
        print(f"Iteration {args.iteration}: Refining HOGs from previous run...")
        
        # Look for the previous run's output
        prev_run_path = f"../HmmRun{args.iteration - 1}"
        # prev_db = f"{prev_run_path}/{args.project_name}HmmDatabase/{args.project_name}fam{args.iteration - 1}"
        prev_db = f"{prev_run_path}/{args.project_name}HmmDatabase/{args.project_name}fam1"
        # super_matrix_path = os.path.join(prev_run_path, base_db_path, 'superMatrix')
        
        # Prepare the directory
        if os.path.exists("./hogDirectory"): 
            shutil.rmtree("./hogDirectory")
            
        # If we have a disagreement matrix/superMatrix, we use it here to filter
        if os.path.exists(prev_run_path + "/hogDirectoryFiltered"):
             print("Copying filtered directory from previous run...")
             shutil.copytree(f"{prev_run_path}/hogDirectoryFiltered", "./hogDirectory")
        else:
             print("Error: Previous filtered directory not found!"); exit(1)

        print("Checking for genes to remove (Disagreements & Dropouts)...")
        genes_to_remove = set()
        
        try:
            sm_path = os.path.join(prev_db, 'superMatrix')
            if os.path.exists(sm_path):
                sm = pd.read_csv(sm_path, sep='\t')
                # Find rows where Hog != HmmHogHit
                bad_rows = sm[sm["Hog"] != sm["HmmHogHit"]]
                genes_to_remove.update(bad_rows["Gene"].tolist())
                print(f"    Found {len(bad_rows)} disagreement genes.")
        except Exception as e:
            print(f"Warning reading superMatrix: {e}")
        
        # Remove genes
        if genes_to_remove:
            clean_genes_from_hogs("./hogDirectory", genes_to_remove)

    # Generate statistics
    print("Generating Statistics (hogMsaStats.tsv)...")
    hog_stats_df = calculate_stats("hogDirectory", "hogMsaStats.tsv")

    if hog_stats_df.empty:
        print("No HOGs found to process.")
        return

    # Filtering logic
    print(f"Filtering HOGs (Min Species >= {args.min_species})...")
    
    output_dir = "hogDirectoryFiltered"
    os.makedirs(output_dir, exist_ok=True)
    
    count_kept = 0
    for index, row in hog_stats_df.iterrows():
        keep_file = False
        # Grass specific filtering
        if args.grass_filter:
            species_count = row['species_count']
            
            # Check presence using the column names which match species abbreviations
            # We use .get(species, 0) to be safe if a species column is missing
            found_pacmad = sum(int(row.get(s, 0)) == 1 for s in pacmad)
            found_bop = sum(int(row.get(s, 0)) == 1 for s in bop)

            if (species_count >= 8) or \
               (found_pacmad == 0 and found_bop >= 6) or \
               (found_bop == 0 and found_pacmad >= 6):
                keep_file = True

        else: 
            # Generic filter: Keep if species count is high enough
            # Can add more complex later
            if row['species_count'] >= args.min_species:
                keep_file = True
        
        if keep_file:
            # Construct filename (handle N0 prefix variations)
            hog_name = str(row['HOG'])
            if not hog_name.startswith("N0."):
                src_file_n0 = os.path.join("hogDirectory", f"N0.{hog_name}.fasta")
                src_file_std = os.path.join("hogDirectory", f"{hog_name}.fasta")
            else:
                src_file_n0 = os.path.join("hogDirectory", f"{hog_name}.fasta")
                src_file_std = src_file_n0

            # Copy file if it exists
            if os.path.exists(src_file_n0):
                shutil.copy(src_file_n0, output_dir)
                count_kept += 1
            elif os.path.exists(src_file_std):
                shutil.copy(src_file_std, output_dir)
                count_kept += 1

    print(f"Filtering Complete. Kept {count_kept} files out of {len(hog_stats_df)}.")
    
    # Recalculate stats for the filtered set
    calculate_stats(output_dir, "hogMsaStats2.tsv")
    print("Done filtering.")

if __name__ == "__main__":
    main()