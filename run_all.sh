#!/bin/bash
#
# Complete workflow script to generate thermodynamic data and calculate equilibrium concentrations
#
# Author: David Lary (davidlary@me.com)
#

set -e  # Exit immediately if a command exits with a non-zero status

echo "============================================================"
echo "Atmospheric Thermodynamic Equilibrium Calculator"
echo "Complete Workflow"
echo "============================================================"

# Check if Docker is available
if command -v docker &> /dev/null; then
    USE_DOCKER=true
    echo "Docker found. Using Docker container for calculations."
else
    USE_DOCKER=false
    echo "Docker not found. Using local Python installation."
    
    # Check if requirements are installed
    if ! pip list | grep -q cantera; then
        echo "Installing required packages..."
        pip install -r requirements.txt
    fi
fi

echo
echo "============================================================"
echo "STAGE 1: Generating thermodynamic data"
echo "============================================================"
echo "Reading species from Species.yaml"
echo "Generating NASA-9 polynomial thermodynamic data..."

if [ "$USE_DOCKER" = true ]; then
    # Build Docker image if it doesn't exist
    if ! docker images | grep -q thermodynamics; then
        echo "Building Docker image..."
        docker build -t thermodynamics .
    fi
    
    # Run using Docker
    docker run -it --rm -v "$(pwd)":/app thermodynamics python thermo_generator.py
else
    # Run locally
    python thermo_generator.py
fi

echo
echo "============================================================"
echo "STAGE 2: Calculating equilibrium concentrations"
echo "============================================================"
echo "Reading configuration from EquilibriumCalculation.yaml"
echo "Using thermodynamic data from Thermodynamics.yaml"
echo "Calculating equilibrium concentrations..."

if [ "$USE_DOCKER" = true ]; then
    # Run using Docker
    docker run -it --rm -v "$(pwd)":/app thermodynamics python EquilibriumCalculation.py
else
    # Run locally
    python EquilibriumCalculation.py
fi

echo
echo "============================================================"
echo "Workflow complete!"
echo "Results are available in the 'results' directory"
echo "============================================================"