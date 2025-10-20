#!/bin/bash

# This script runs the calculate_dielectric_derivatives.py script over a user-defined
# range of laser energies, automating the input prompt for each run.
# It will generate a separate set of output files (e.g., dielectric_derivatives_1.00)
# for each energy step.

# --- Get User Input ---
echo "Enter the start energy in eV (e.g., 1.0):"
read START_EV

echo "Enter the end energy in eV (e.g., 5.0):"
read END_EV

echo "Enter the energy increment in eV (e.g., 0.1):"
read INCREMENT_EV

echo ""
echo "Starting energy scan from ${START_EV} eV to ${END_EV} eV in steps of ${INCREMENT_EV} eV."
echo "----------------------------------------------------"

# --- Main Loop ---
# We use the 'seq' command to generate the sequence of floating-point numbers.
for energy in $(seq ${START_EV} ${INCREMENT_EV} ${END_EV})
do
  echo ""
  echo "--- Running analysis for ${energy} eV ---"
  
  # The 'echo' command sends the current energy value into the standard
  # input of the Python script, answering the prompt automatically.
  echo "${energy}" | python calculate_dielectric_derivatives.py
  
  # Optional: Check if the python script ran successfully
  if [ $? -ne 0 ]; then
    echo "Error running Python script for ${energy} eV. Aborting scan."
    exit 1
  fi
done

echo ""
echo "----------------------------------------------------"
echo "Energy scan complete."