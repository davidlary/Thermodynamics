#!/bin/bash
set -e

# Check if the first argument is run-all
if [ "$1" = "run-all" ]; then
    echo "========================================"
    echo "Running complete thermodynamic workflow"
    echo "========================================"
    
    # Step 1: Generate thermodynamic data
    echo -e "\n[Step 1] Generating thermodynamic data..."
    python thermo_generator.py
    
    # Step 2: Calculate equilibrium concentrations
    echo -e "\n[Step 2] Calculating equilibrium concentrations..."
    python EquilibriumCalculation.py
    
    echo -e "\nWorkflow completed successfully!"
    
    # List results
    echo -e "\nResults:"
    ls -l results/
    
    # Show log files
    echo -e "\nLog files:"
    ls -l logs/
    
# If the command starts with "python", run Python with the arguments
elif [ "${1:0:6}" = "python" ]; then
    exec "$@"
    
# Otherwise pass the arguments to bash
else
    exec "$@"
fi