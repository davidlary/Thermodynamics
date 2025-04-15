#!/bin/bash
#
# Complete workflow script to generate thermodynamic data and calculate equilibrium concentrations
#
# Author: David Lary (davidlary@me.com)
#

set -e  # Exit immediately if a command exits with a non-zero status

# Function to print a section header
print_header() {
    local title="$1"
    echo
    echo "============================================================"
    echo "$title"
    echo "============================================================"
}

# Main header
print_header "Atmospheric Thermodynamic Equilibrium Calculator - Complete Workflow"

# Check if logs directory exists, create if not
if [ ! -d "logs" ]; then
    mkdir -p logs
    echo "Created logs directory."
fi

# Check if results directory exists, create if not
if [ ! -d "results" ]; then
    mkdir -p results/plots results/summary
    echo "Created results directory structure."
fi

# Check if Docker is available
if command -v docker &> /dev/null; then
    USE_DOCKER=true
    echo "Docker found. Using Docker container for calculations."
    
    # Check if thermodynamics image exists, build if not
    if ! docker images | grep -q thermodynamics; then
        echo "Building Docker image..."
        docker build -t thermodynamics .
    fi
    
    # Run the complete workflow in Docker with a single command
    print_header "Running complete workflow in Docker container"
    docker run -it --rm -v "$(pwd)":/app thermodynamics run-all
    
    # Check if the execution was successful
    if [ $? -eq 0 ]; then
        print_header "Workflow completed successfully in Docker container!"
        echo "Results are available in the 'results' directory"
        echo "Logs are available in the 'logs' directory"
    else
        print_header "ERROR: Workflow failed in Docker container"
        echo "Check logs in the 'logs' directory for details"
        exit 1
    fi
    
else
    USE_DOCKER=false
    echo "Docker not found. Using local Python installation."
    
    # Check if requirements are installed
    if ! pip list | grep -q cantera; then
        echo "Installing required packages..."
        pip install -r requirements.txt
    fi
    
    # Run Stage 1: Generate thermodynamic data
    print_header "STAGE 1: Generating thermodynamic data"
    echo "Reading species from Species.yaml"
    echo "Generating NASA-9 polynomial thermodynamic data..."
    python thermo_generator.py
    
    # Run Stage 2: Calculate equilibrium concentrations
    print_header "STAGE 2: Calculating equilibrium concentrations"
    echo "Reading configuration from EquilibriumCalculation.yaml"
    echo "Using thermodynamic data from Thermodynamics.yaml"
    echo "Calculating equilibrium concentrations..."
    python EquilibriumCalculation.py
    
    # Check if the execution was successful
    if [ $? -eq 0 ]; then
        print_header "Workflow completed successfully!"
        echo "Results are available in the 'results' directory"
        echo "Logs are available in the 'logs' directory"
    else
        print_header "ERROR: Workflow failed"
        echo "Check logs in the 'logs' directory for details"
        exit 1
    fi
fi

# Display logs and results
echo
echo "Latest log files:"
ls -lt logs/ | head -5

echo
echo "Results:"
ls -lt results/ | head -5