#!/bin/bash
ABSOLUTE_PATH_OF_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Specify the directory path
directory="$ABSOLUTE_PATH_OF_SCRIPT"/../

# Use find to locate files without an extension in subdirectories, excluding .git directory
find "$directory" -type f ! -path "$directory/.git/*" ! -name "*.*" -exec rm -f {} \;


echo "Executable files removed successfully."

