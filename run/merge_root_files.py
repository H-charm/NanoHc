#!/usr/bin/env python3
import os
import subprocess
import sys

def merge_root_files(dir1, dir2, output_dir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # List all ROOT files in the first directory
    files1 = [f for f in os.listdir(dir1) if f.endswith('.root')]
    
    # Loop over each file in dir1, check for its counterpart in dir2, then merge
    for filename in files1:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if os.path.exists(file2):
            output_file = os.path.join(output_dir, filename)
            # Build the haddnano.py command: output file then the two input files
            cmd = ["haddnano.py", output_file, file1, file2]
            print(f"Merging {filename} ...")
            result = subprocess.run(cmd)
            if result.returncode != 0:
                print(f"Error merging {filename}.")
        else:
            print(f"Warning: {filename} not found in {dir2}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 merge_root_files.py <dir1> <dir2> <output_dir>")
        sys.exit(1)
    
    directory1 = sys.argv[1]
    directory2 = sys.argv[2]
    output_directory = sys.argv[3]
    
    merge_root_files(directory1, directory2, output_directory)
