import os
import sys

def check_if_dir_exists(path):
    if os.path.exists(path):
        while True:  # Loop until valid input is received
            user_input = input(f"The directory '{path}' already exists. Do you want to overwrite it? (y/n): ").lower()
            if user_input in ('y', 'n'):
                break  # Exit the loop if input is valid
            print("Invalid input. Please enter 'y' or 'n'.")
        if user_input == 'n': sys.exit(0)
        
## to divide datasets by number of files per job
def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]