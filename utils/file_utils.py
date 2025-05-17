"""
DOLFINx Modular Structural Solver - File utilities

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
def load_config(file_path):
    import yaml
    if file_path.endswith('.yaml') or file_path.endswith('.yml'):
        with open(file_path) as f:
            return yaml.safe_load(f)
    else:
        raise ValueError("Unsupported file type. Please use .yaml or .yml files.")
    
if __name__ == '__main__':
    print("This module is not meant to be run directly.")