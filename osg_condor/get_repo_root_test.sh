#!/bin/bash

#chatgpt assist with method to derive the path to the github repository for use with shell scripts to call other scripts in this repository

# Save the current working directory
original_dir=$(pwd)
# Get the directory where the script is located
script_dir=$(dirname "$(realpath "$0")")
# Change to the script's directory
cd "$script_dir"
# Find the root of the GitHub repository
repo_root=$(git rev-parse --show-toplevel)
# Return to the original directory
cd "$original_dir"
# Now `repo_root` holds the path to the root of the GitHub repository
echo "The root of the GitHub repository is: $repo_root"
