#!/bin/bash
ABSOLUTE_PATH_OF_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Specify the directory path
directory="$ABSOLUTE_PATH_OF_SCRIPT"/../

# Use find to locate files without an extension in subdirectories, excluding .git directory
find "$directory" -type f ! -path "$directory/.git/*" ! -name "*.*" -exec rm -f {} \;

# Use find to locate and remove directories with a specific extension (.dSYM)
find "$directory" -type d -name "*.dSYM" -exec rm -rf {} \;

echo "Executable files removed successfully."

