"""
DOLFINx Modular Structural Solver - Main entry point

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
import argparse
import os

from src.solver import run_modal_analysis

# get the config path from the command line
def parse_args():
    parser = argparse.ArgumentParser(description="Run modal analysis with FEniCS.")
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to the configuration file (YAML format).",
    )
    args = parser.parse_args()
    if not os.path.exists(args.config):
        parser.error(f"Configuration file '{args.config}' does not exist.")
    return args

def main():
    args = parse_args()
    config_path = args.config
    run_modal_analysis(config_path)

if __name__ == "__main__":
    main()