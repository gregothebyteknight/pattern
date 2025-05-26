#!/bin/bash

# filepath: ~/Documents/shapes/code/multiply.sh

# Base directory containing the simulation folders
BASE_DIR="../../spheres/data"

# Loop through simulation numbers 1 to 30
for i in $(seq 1 30); do
  # Format the folder name
  FOLDER_NAME="sim_$i"
  SIM_FOLDER="$BASE_DIR/$FOLDER_NAME"

  # Skip if the folder does not exist (note spaces around [ and !)
  if [ ! -d "$SIM_FOLDER" ]; then
    echo "Folder $SIM_FOLDER does not exist. Skipping..."
    continue
  fi

  # Find the input file (assuming one CSV file per folder)
  INPUT_FILE=$(find "$SIM_FOLDER" -type f -name "*.csv" | head -n 1)

  # Skip if no CSV found
  if [ ! -f "$INPUT_FILE" ]; then
    echo "No CSV file found in $SIM_FOLDER. Skipping..."
    continue
  fi

  echo "Processing folder: $SIM_FOLDER"
  echo "Running angle.r for input file: $INPUT_FILE"

  # Run the angle.r script and capture its output
  ANGLE_OUT=$(Rscript ~/Documents/shapes/code/angle.r "$INPUT_FILE")
  echo "$ANGLE_OUT"

  # Extract Optimal max_r from angle.r output
  RAW_MAX=$(printf '%s\n' "$ANGLE_OUT" \
            | grep -i "Optimal max_r")

  # extract the numeric portion (handles decimal)  
  NUM=$(printf '%s\n' "$RAW_MAX" \
        | sed -E 's/.*: *([0-9]+(\.[0-9]+)?).*/\1/')

  # multiply by 1.1 
  MAX_R=$(echo "$NUM * 1.1" | bc -l) 

  # Skip slice.r if parsing failed
  if [ -z "$MAX_R" ]; then
    echo "  → Could not parse Optimal max_r; skipping slice.r"
    continue
  fi

  echo " → Parsed raw max_r = $NUM; scaled ×1.1 → $MAX_R"
  echo " → Running slice.r with r_max = $MAX_R"
  Rscript ~/Documents/shapes/code/slice.r "$INPUT_FILE" "$MAX_R"
done
