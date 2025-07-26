import sys, os
import argparse
from edlib import align
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from tqdm import tqdm
import tempfile

def read_fasta(fname):
    In = open(fname)
    temp = In.read().split("\n")
    In.close()
    fasta = {}
    for row in temp:
        if row.find(">")!=-1:
            active_seq = row[1:]
            fasta[active_seq] = ""
        else:
            fasta[active_seq] += row
    return fasta

def separate_fasta(fasta):
    ID_list = [ID for ID in fasta]
    seqlist = [fasta[ID] for ID in ID_list]
    return ID_list, seqlist

def is_fastq(file_path):
    """Check if a file is in FASTQ format."""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('@'):
                return False
            
            # Read next 2 lines
            _ = f.readline()  # Sequence
            third_line = f.readline().strip()
            
            # Check if third line starts with '+'
            return third_line.startswith('+')
    except:
        return False

def convert_fastq_to_fasta(file_path):
    """Convert a FASTQ file to FASTA format, discarding quality scores."""
    if not is_fastq(file_path):
        return file_path
    
    temp_fasta = file_path + ".temp.fasta"
    try:
        with open(file_path, 'r') as f_in, open(temp_fasta, 'w') as f_out:
            while True:
                header = f_in.readline()
                if not header:
                    break
                
                sequence = f_in.readline()
                _ = f_in.readline()  # '+' line
                _ = f_in.readline()  # Quality scores
                
                if header and sequence:
                    f_out.write('>' + header[1:])
                    f_out.write(sequence)
        return temp_fasta
    except Exception as e:
        print(f"Error converting FASTQ to FASTA: {e}")
        # If conversion fails, return the original file path
        return file_path

def helper_single_file(i):
    """Helper function for pairwise comparisons within a single file."""
    global ID_list, seqlist, temp_dir, mode, query_is_first
    
    temp_out = temp_dir + str(i) + "_" + str(os.getpid()) + ".ijk"
    Out = open(temp_out, 'w+')
    for j in range(i+1, len(ID_list)):
        # CRITICAL: Order matters for HW mode!
        if mode == "HW":
            if query_is_first:
                # sequence i is the query, sequence j is the target
                # query (i) can be infix-ed into target (j)
                result = align(seqlist[i], seqlist[j], mode=mode)
            else:
                # sequence j is the query, sequence i is the target
                # query (j) can be infix-ed into target (i)
                result = align(seqlist[j], seqlist[i], mode=mode)
        else:
            # In NW mode: global alignment, order doesn't matter
            result = align(seqlist[i], seqlist[j], mode=mode)
            
        Out.write(ID_list[i] + "\t" + ID_list[j] + "\t" + str(result['editDistance']) + "\n")
    Out.close()
    return temp_out

def helper_two_files(args):
    """Helper function for comparisons between two different files.
       For each sequence in chunk from input2, compare to all sequences in input1."""
    chunk_idx, read_ids, read_seqs = args
    global ID_list1, seqlist1, temp_dir, mode, query_is_first, find_closest
    
    # Create two output files: one for all distances, one for tracking closest matches
    temp_out = temp_dir + "all_" + str(chunk_idx) + "_" + str(os.getpid()) + ".ijk"
    
    # If we need to find closest matches
    if find_closest:
        temp_closest = temp_dir + "closest_" + str(chunk_idx) + "_" + str(os.getpid()) + ".ijk"
        closest_out = open(temp_closest, 'w+')
    
    with open(temp_out, 'w+') as out_file:
        for i, (read_id, read_seq) in enumerate(zip(read_ids, read_seqs)):
            # For each read, find distances to all references
            min_dist = float('inf')
            min_refs = []
            
            for j, ref_seq in enumerate(seqlist1):
                # CRITICAL: Order matters for HW mode!
                if mode == "HW":
                    if query_is_first:
                        # ref is the query, read is the target
                        # query (ref) can be infix-ed into target (read)
                        result = align(ref_seq, read_seq, mode=mode)
                    else:
                        # read is the query, ref is the target
                        # query (read) can be infix-ed into target (ref)
                        result = align(read_seq, ref_seq, mode=mode)
                else:
                    # In NW mode: global alignment, order doesn't matter
                    result = align(ref_seq, read_seq, mode=mode)
                
                dist = result['editDistance']
                out_file.write(ID_list1[j] + "\t" + read_id + "\t" + str(dist) + "\n")
                
                # Track minimum distance for closest match output
                if find_closest:
                    if dist < min_dist:
                        min_dist = dist
                        min_refs = [ID_list1[j]]
                    elif dist == min_dist:
                        min_refs.append(ID_list1[j])
            
            # Write closest match information
            if find_closest:
                closest_out.write(read_id + "\t" + str(min_dist) + "\t" + "\t".join(min_refs) + "\n")
    
    if find_closest:
        closest_out.close()
        return (temp_out, temp_closest)
    else:
        return temp_out

def print_alignment_details():
    """Print detailed explanation of alignment modes."""
    print("\nAlignment Mode Details:")
    print("----------------------")
    print("NW (Global): Standard edit distance calculation. Both sequences are aligned from end to end.")
    print("              Every mismatch, insertion, or deletion counts towards the edit distance.")
    print()
    print("HW (Infix):  Query sequence can align to any part of the target sequence.")
    print("              Gaps at the beginning and end of the QUERY are NOT penalized.")
    print("              This makes a huge difference in the result! For example:")
    print("                align('TTT', 'AAATTTAAA', mode='HW') → edit distance = 0")
    print("                align('AAATTTAAA', 'TTT', mode='HW') → edit distance = 6")
    print()
    print("When using --query_is_first=True (default for HW mode):")
    print("  - For two-file mode:     Sequences from file1 are QUERIES that get infix-ed into sequences from file2")
    print("  - For single-file mode:  Each sequence i is the QUERY that gets infix-ed into each sequence j")
    print()
    print("When using --query_is_first=False:")
    print("  - For two-file mode:     Sequences from file2 are QUERIES that get infix-ed into sequences from file1")
    print("  - For single-file mode:  Each sequence j is the QUERY that gets infix-ed into each sequence i")
    print("----------------------")

def split_into_chunks(id_list, seq_list, chunk_size):
    """Split lists into chunks for parallel processing."""
    for i in range(0, len(id_list), chunk_size):
        yield (i // chunk_size, id_list[i:i + chunk_size], seq_list[i:i + chunk_size])

def combine_files(file_list, output_file, file_type="distance"):
    """Combine multiple files into one output file with error handling."""
    print(f"Combining {len(file_list)} {file_type} files...")
    
    successful_files = 0
    missing_files = []
    total_lines = 0
    
    try:
        with open(output_file, 'w') as fout:
            for fname in file_list:
                if not os.path.exists(fname):
                    missing_files.append(fname)
                    print(f"[WARNING] Missing expected temp file: {fname}")
                    continue
                
                try:
                    with open(fname, 'r') as fin:
                        lines_in_file = 0
                        for line in fin:
                            fout.write(line)
                            lines_in_file += 1
                        total_lines += lines_in_file
                        successful_files += 1
                except Exception as e:
                    print(f"[ERROR] Failed to read file {fname}: {e}")
                    continue
        
        print(f"Successfully combined {successful_files}/{len(file_list)} files")
        print(f"Total lines written: {total_lines}")
        
        if missing_files:
            print(f"[WARNING] {len(missing_files)} files were missing")
            if len(missing_files) <= 10:
                for mf in missing_files:
                    print(f"  - {mf}")
            else:
                print(f"  - First 10 missing files:")
                for mf in missing_files[:10]:
                    print(f"    - {mf}")
                print(f"  - ... and {len(missing_files) - 10} more")
        
        # Verify output file was created
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"Output file created: {output_file} (size: {file_size} bytes)")
        else:
            print(f"[ERROR] Output file was not created: {output_file}")
            
    except Exception as e:
        print(f"[ERROR] Failed to create output file {output_file}: {e}")
        raise

def safe_remove_files(file_list):
    """Safely remove temporary files with error handling."""
    print(f"Cleaning up {len(file_list)} temporary files...")
    
    removed = 0
    failed = 0
    
    for fname in file_list:
        try:
            if os.path.exists(fname):
                os.remove(fname)
                removed += 1
            else:
                # File already doesn't exist, which is fine
                pass
        except Exception as e:
            print(f"[WARNING] Failed to remove temp file {fname}: {e}")
            failed += 1
    
    print(f"Removed {removed} files, {failed} failures")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate pairwise edit distances between sequences using parallel processing.')
    
    # Required parameters
    parser.add_argument('--input1', required=True, help='First input file (FASTA or FASTQ format)')
    parser.add_argument('--output', required=True, help='Output file for all pairwise distances (.ijk format)')
    parser.add_argument('--temp_dir', required=True, help='Directory for temporary files')
    
    # Optional parameters
    parser.add_argument('--input2', help='Second input file (optional, for cross-comparisons)')
    parser.add_argument('--threads', type=int, default=-1, help='Number of threads (-1 for automatic)')
    parser.add_argument('--mode', choices=['NW', 'HW'], default='NW', 
                        help='Alignment mode: NW (global) or HW (infix)')
    parser.add_argument('--query_is_first', type=bool, default=True,
                        help='For HW mode: if True, the first sequence is the query (gaps at its ends not penalized)')
    parser.add_argument('--closest_output', help='Output file for closest matches. When specified, '
                        'for each sequence in input2, find the closest sequence(s) in input1')
    parser.add_argument('--chunk_size', type=int, default=100, 
                        help='Number of sequences to process in each parallel chunk (default: 100)')
    parser.add_argument('--details', action='store_true', help='Show detailed explanation of alignment modes')
    
    args = parser.parse_args()
    
    if args.details:
        print_alignment_details()
        sys.exit(0)
    
    # Set global variables
    global ID_list, seqlist, ID_list1, seqlist1, ID_list2, seqlist2, temp_dir, mode, query_is_first, find_closest
    
    # Set variables from arguments
    temp_dir = args.temp_dir
    mode = args.mode
    query_is_first = args.query_is_first
    find_closest = args.closest_output is not None
    
    # Ensure temp_dir ends with a slash
    if temp_dir[-1] != '/':
        temp_dir += '/'
    
    # Create temp directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)
    
    # Handle thread count
    n_threads = args.threads
    if n_threads < 1 or n_threads > cpu_count()-1:
        n_threads = max(1, cpu_count()-1)
    
    # Pre-process input files if they're FASTQ
    in_fname1 = convert_fastq_to_fasta(args.input1)
    temp_files = []
    if in_fname1 != args.input1:
        temp_files.append(in_fname1)
    
    # Print alignment mode info
    print(f"Using alignment mode: {mode}")
    if mode == "HW":
        if query_is_first:
            print("HW mode with query_is_first=True: gaps at the ends of the FIRST sequence are NOT penalized")
        else:
            print("HW mode with query_is_first=False: gaps at the ends of the SECOND sequence are NOT penalized")
    
    if args.input2:
        # Two-file mode
        in_fname2 = convert_fastq_to_fasta(args.input2)
        if in_fname2 != args.input2:
            temp_files.append(in_fname2)
        
        # Process both files
        ID_list1, seqlist1 = separate_fasta(read_fasta(in_fname1))
        ID_list2, seqlist2 = separate_fasta(read_fasta(in_fname2))
        
        print(f"Processing {len(ID_list1)} sequences from file 1 against {len(ID_list2)} sequences from file 2")
        print(f"Using {n_threads} threads")
        
        if mode == "HW":
            if query_is_first:
                print("Sequences from file 1 will be infix-ed (as queries) into sequences from file 2")
                if find_closest:
                    print(f"Finding closest amplicon(s) from file 1 for each read in file 2")
                    print(f"Results will be written to {args.closest_output}")
            else:
                print("Sequences from file 2 will be infix-ed (as queries) into sequences from file 1")
                if find_closest:
                    print(f"Finding closest amplicon(s) from file 1 for each read in file 2")
                    print(f"Results will be written to {args.closest_output}")
        
        # Split input2 (reads) into chunks for parallel processing
        chunks = list(split_into_chunks(ID_list2, seqlist2, args.chunk_size))
        total_chunks = len(chunks)
        
        print(f"Processing data in {total_chunks} chunks")
        
        # Process in parallel with a progress bar
        all_fnames = []
        closest_fnames = []
        
        with Pool(processes=n_threads) as pool:
            results = list(tqdm(pool.imap(helper_two_files, chunks), total=total_chunks, desc="Processing chunks"))
            
            # Separate all distances and closest matches files
            if find_closest:
                for result in results:
                    if result is not None:
                        if isinstance(result, tuple):
                            all_fname, closest_fname = result
                            all_fnames.append(all_fname)
                            closest_fnames.append(closest_fname)
                        else:
                            all_fnames.append(result)
            else:
                all_fnames = [r for r in results if r is not None]
        
        # Combine all distances results
        combine_files(all_fnames, args.output, "distance")
        
        # Combine closest matches results if needed
        if find_closest and closest_fnames:
            combine_files(closest_fnames, args.closest_output, "closest match")
        
        # Clean up temporary files
        print("Cleaning up temporary files...")
        all_temp_files = all_fnames[:]
        if find_closest:
            all_temp_files.extend(closest_fnames)
        
        safe_remove_files(all_temp_files)
        
    else:
        # Single-file mode
        ID_list, seqlist = separate_fasta(read_fasta(in_fname1))
        
        print(f"Processing {len(ID_list)} sequences for all pairwise comparisons")
        print(f"Using {n_threads} threads")
        
        if mode == "HW":
            if query_is_first:
                print("Each sequence i will be infix-ed (as query) into each sequence j (where j > i)")
            else:
                print("Each sequence j will be infix-ed (as query) into each sequence i (where j > i)")
        
        # Create a list of indices to process
        indices = list(range(len(ID_list)-1))
        
        # Process in parallel with a progress bar
        with Pool(processes=n_threads) as pool:
            fnames = list(tqdm(pool.imap(helper_single_file, indices), total=len(indices), desc="Processing pairs"))
        
        # Filter out None results
        fnames = [f for f in fnames if f is not None]
        
        # Combine results
        combine_files(fnames, args.output, "distance")
        
        # Clean up temporary files
        print("Cleaning up temporary files...")
        safe_remove_files(fnames)
    
    # Clean up converted files
    for temp_file in temp_files:
        try:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        except Exception as e:
            print(f"[WARNING] Failed to remove converted temp file {temp_file}: {e}")
    
    # Final verification
    if os.path.exists(args.output):
        print(f"✓ Done! Results written to {args.output}")
    else:
        print(f"✗ ERROR: Output file was not created: {args.output}")
        
    if find_closest:
        if args.closest_output and os.path.exists(args.closest_output):
            print(f"✓ Closest matches written to {args.closest_output}")
        else:
            print(f"✗ ERROR: Closest matches file was not created: {args.closest_output}")
