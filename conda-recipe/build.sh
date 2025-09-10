# conda-recipe/build.sh
#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Make sure that the bin directory exists
mkdir -p $PREFIX/bin

# Copy the processing script to the bin directory
cp ../scripts/process_cram_files.sh $PREFIX/bin/

# Make the script executable
chmod +x $PREFIX/bin/process_cram_files.sh

# If there are other executables, copy them similarly

