FROM python:3.10-slim

LABEL maintainer="David Lary <davidlary@me.com>"
LABEL description="Atmospheric Thermodynamic Equilibrium Calculator"

# Install system dependencies required for Cantera
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    swig \
    libboost-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create and set working directory
WORKDIR /app

# Copy requirements file
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY . .

# Create directories for outputs
RUN mkdir -p results/plots results/summary

# Set environment variables
ENV PYTHONUNBUFFERED=1

# Command to run when container starts
CMD ["bash"]

# Usage instructions as comment
# To generate thermodynamic data:
#   docker run -it --rm -v $(pwd):/app thermodynamics python thermo_generator.py
# To calculate equilibrium concentrations:
#   docker run -it --rm -v $(pwd):/app thermodynamics python EquilibriumCalculation.py