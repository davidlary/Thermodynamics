# Atmospheric Thermodynamic Equilibrium Calculator

This project calculates thermodynamic equilibrium concentrations of chemical species as a function of temperature and pressure using the Cantera library, with a special focus on accuracy in the cold temperature range (100-400K) for atmospheric research applications.

Developed by David Lary (davidlary@me.com)

## Process Overview

The calculation is performed in two distinct stages:

### Stage 1: Generate Thermodynamic Data
First, we generate high-quality NASA-9 polynomial thermodynamic data from various sources:
```bash
# Standard method
python thermo_generator.py

# Using Docker
docker run -it --rm -v $(pwd):/app thermodynamics python thermo_generator.py
```
This reads the list of species from `Species.yaml` and generates `Thermodynamics.yaml` with NASA-9 polynomial data.

### Stage 2: Calculate Equilibrium Concentrations
Then, we use the generated thermodynamic data to calculate equilibrium concentrations:
```bash
# Standard method
python EquilibriumCalculation.py

# Using Docker
docker run -it --rm -v $(pwd):/app thermodynamics python EquilibriumCalculation.py
```
This reads configuration from `EquilibriumCalculation.yaml` and uses the thermodynamic data from `Thermodynamics.yaml`.

IMPORTANT: The two stages must be run in sequence. You must first generate the thermodynamic data (Stage 1) before calculating equilibrium concentrations (Stage 2).

### Automated Workflow

For convenience, we provide a shell script that runs both stages in sequence:
```bash
./run_all.sh
```

This script automatically detects if Docker is available and uses it if possible, otherwise falls back to the local Python installation.

## Key Components
1. **Thermodynamic Data Generator** (`thermo_generator.py`): Creates NASA-9 polynomial thermodynamic data from multiple sources, with special focus on cold temperature accuracy for atmospheric research
2. **Equilibrium Concentration Calculator** (`EquilibriumCalculation.py`): Calculates equilibrium concentrations across a temperature range using Cantera, with our custom thermodynamic data

## Configuration Files
- `Species.yaml`: List of chemical species to include in the thermodynamic data generation
- `Thermodynamics.yaml`: Generated thermodynamic data in NASA-9 polynomial format
- `EquilibriumCalculation.yaml`: Configuration for equilibrium calculations (temperature range, pressure, initial mixture, etc.)

## Features

- Reads chemical species from a YAML file
- Generates NASA-9 polynomial thermodynamic data for Cantera
- Prioritized thermodynamic data lookup: Burcat → CEA → NASA → NIST
- Theoretical first principles calculations when no experimental data exists
- Local caching of thermodynamic data for improved performance
- Enhanced temperature range support from 100K to 6000K
- Optimized accuracy in cold temperature range (100-400K) for atmospheric research
- Robust temperature range transition continuity verification and optimization
- Support for all relevant chemical elements and charged species
- No placeholder data - only real or theoretically calculated data from validated sources

## Requirements

- Python 3.6+
- Cantera
- PyYAML
- Requests

## Installation

### Method 1: Standard Installation

```bash
pip install -r requirements.txt
```

### Method 2: Docker (Recommended for Reproducibility)

We provide a Docker container that comes pre-loaded with all dependencies:

1. Build the Docker image:
```bash
docker build -t thermodynamics .
```

2. Run the container with your local directory mounted:
```bash
docker run -it --rm -v $(pwd):/app thermodynamics
```

## Usage

### 1. Define Chemical Species

Define your species in `Species.yaml`:

```yaml
# List of chemical species
species:
  - N2
  - O2
  - H2O
  - CO2
  # Add more species as needed
```

### 2. Generate Thermodynamic Data

Generate NASA-9 polynomial thermodynamic data:

```bash
# Standard method
python thermo_generator.py

# Using Docker
docker run -it --rm -v $(pwd):/app thermodynamics python thermo_generator.py
```

This will create a `Thermodynamics.yaml` file with the NASA-9 polynomial data ready for use in Cantera.

### 3. Configure Equilibrium Calculations

Customize your equilibrium calculation parameters in `EquilibriumCalculation.yaml`:

```yaml
# Input Data
input:
  thermo_data_file: "Thermodynamics.yaml"  # Path to thermodynamic data file
  use_cantera_database: false  # Use our custom thermodynamic data from Thermodynamics.yaml

# Calculation Parameters
calculation:
  temperature_range: [200, 2000]  # K (min, max)
  log_temperature: true  # Use logarithmic spacing for better low-temp resolution
  temperature_points: 100  # Number of temperature points to calculate
  pressure: 101325  # Pa (1 atm)

# Initial Mixture Composition (mole fractions)
initial_mixture:
  N2: 0.78084
  O2: 0.20946
  Ar: 0.00934
  CO2: 0.00041
  H2O: 0.00500  # ~50% relative humidity
  # Add more species with their mole fractions
```

IMPORTANT NOTE: The equilibrium calculator is designed to use our custom thermodynamic data from Thermodynamics.yaml, which provides superior accuracy especially in the cold temperature range (100-400K). You must run thermo_generator.py first to generate this file before running the equilibrium calculator.

### 4. Run Equilibrium Calculations

Calculate equilibrium concentrations across the temperature range:

```bash
# Standard method
python EquilibriumCalculation.py

# Using Docker
docker run -it --rm -v $(pwd):/app thermodynamics python EquilibriumCalculation.py
```

This will:
- Calculate equilibrium concentrations at each temperature point
- Generate concentration vs. temperature plots for each species
- Create summary plots with multiple species per page
- Save all results to CSV files for further analysis

## How It Works

The system follows these steps:

1. Reads the list of species from `Species.yaml`
2. Extracts all unique elements from the chemical formulas
3. Downloads and caches the complete Burcat thermodynamic database
4. For each species, retrieves thermodynamic data in this priority order:
   - Burcat's Thermodynamic Database (cached locally, updated weekly)
   - NASA CEA Database
   - NASA Database
   - NIST Database
   - Theoretical calculations from first principles (if no experimental data exists)
5. Locally caches all retrieved data to speed up future runs
6. Performs temperature range continuity checks and optimizations
7. Applies targeted optimizations for the cold temperature range (100-400K)
8. Generates a Cantera-compatible YAML file with the complete thermodynamic data
9. Includes detailed source attribution and validation metadata for each species
10. Creates a data sources report summarizing which database was used for each species

### Thermodynamic Data Sources

The system prioritizes data sources in this order:

1. **Burcat Database**: The most comprehensive and up-to-date source of thermodynamic data, maintained by Alexander Burcat at the Technion - Israel Institute of Technology.
2. **NASA CEA Database**: Chemical Equilibrium with Applications database from NASA.
3. **NASA Thermodynamic Database**: General NASA thermodynamic data.
4. **NIST Chemistry WebBook**: National Institute of Standards and Technology data.
5. **Theoretical First Principles**: When no experimental data exists, uses statistical thermodynamics and molecular properties to calculate thermodynamic data theoretically.

The entire Burcat database is downloaded and cached locally, with automatic version checking to ensure the data remains current without unnecessary downloads.

### Theoretical First Principles Calculations

For species without experimental data, the system uses rigorous theoretical approaches:

- **Statistical Thermodynamics**: Calculates properties based on molecular structure
- **Degrees of Freedom Analysis**: Accounts for translational, rotational, and vibrational modes
- **Molecular Structure Recognition**: Identifies linear molecules, monatomic species, and charged ions
- **Temperature-Dependent Properties**: Generates appropriate polynomial coefficients for all temperature ranges
- **Validation Metadata**: Clearly identifies theoretically derived data with confidence levels

This approach ensures all species have scientifically sound thermodynamic data without resorting to arbitrary placeholders.

### NASA-9 Polynomial Format

This project uses the NASA-9 polynomial format for thermodynamic data which provides:

- Multiple temperature ranges for each species:
  - Cold range (100-400K) for atmospheric research applications
  - Low range (typically 300-1000K)
  - High range (typically 1000-6000K)
- Nine coefficients (a₁-a₉) for each temperature range
- These coefficients define polynomial expressions for:
  - Heat capacity: Cp°/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴
  - Enthalpy: H°/RT = a₁ + (a₂/2)T + (a₃/3)T² + (a₄/4)T³ + (a₅/5)T⁴ + a₆/T
  - Entropy: S°/R = a₁ln(T) + a₂T + (a₃/2)T² + (a₄/3)T³ + (a₅/4)T⁴ + a₇

Where T is temperature in Kelvin and R is the gas constant.

### Temperature Range Optimization

The system employs sophisticated optimization techniques to ensure thermodynamic accuracy across all temperature ranges:

1. **Cold Range Optimization (100-400K)**: Specially optimized coefficients for atmospheric research applications, using iterative refinement of thermodynamic properties with higher weight given to lower temperatures.

2. **Transition Continuity Verification**: Checks for discontinuities in Cp, H, and S at all temperature range transition points.

3. **Transition Continuity Optimization**: Automatically adjusts coefficients to minimize discontinuities while preserving real thermodynamic data accuracy.

4. **Validation Metadata**: All output includes detailed validation information including transition diagnostics and optimization details.

## Equilibrium Concentration Calculator

The `EquilibriumCalculation.py` script provides powerful capabilities for calculating chemical equilibrium concentrations using Cantera and the NASA-9 thermodynamic data:

### Features

- **Temperature Range Scanning**: Calculates equilibrium at each point across a specified temperature range
- **Customizable Initial Mixture**: Define arbitrary initial gas compositions with precise mole fractions
- **Comprehensive Output**: Generates both numerical data (CSV) and visual representations (plots)
- **Species Filtering**: Focus on specific species of interest or exclude others
- **Visualization Options**: Individual species plots and multi-species summary plots
- **Performance Optimization**: Logarithmic temperature distribution for better resolution at lower temperatures
- **Detailed Logging**: Tracks progress, initial conditions, and maximum concentrations
- **Custom Thermodynamic Data**: Uses the high-quality NASA-9 polynomial data from our thermodynamic data generator
- **Cold Range Accuracy**: Benefits from the specialized cold temperature range optimizations (100-400K)
- **Charged Species Support**: Full support for ions and charged species in equilibrium calculations
- **Robust Error Handling**: Gracefully handles missing species and calculation failures
- **Format Conversion**: Automatically converts our NASA-9 format to Cantera-compatible format

### Configuration Options

The `EquilibriumCalculation.yaml` file provides extensive configuration options:

- **Input Data**:
  - Uses our custom NASA-9 polynomial thermodynamic data from Thermodynamics.yaml
  - Provides superior accuracy especially in the 100-400K range for atmospheric applications
- **Temperature Range**:
  - Define minimum and maximum temperatures and number of points
  - Option for logarithmic or linear temperature spacing
- **Pressure**: Set pressure for equilibrium calculations
- **Solver Parameters**: Customize tolerance and maximum iterations for the equilibrium solver
- **Output Options**:
  - Control file formats, plotting styles, and what information to log
  - Specify plot scales (linear or logarithmic)
  - By default, plots all species (can be filtered using focus_species if needed)
- **Initial Mixture**: Define the starting gas composition with precise mole fractions

### Output Files

The calculator generates several types of output:

1. **CSV Data File**: Contains all equilibrium concentrations at each temperature
2. **Individual Species Plots**: Concentration vs. temperature for each species
3. **Summary Plots**: Multiple species per page, organized by concentration
4. **Log File**: Records calculation progress and important results

## Applications

This tool is particularly useful for:

- **Atmospheric Research**: Enhanced accuracy in the 100-400K range makes it ideal for stratospheric and mesospheric chemistry studies
- **Chemical Engineering**: Accurate prediction of equilibrium concentrations across wide temperature ranges
- **Combustion Analysis**: High-temperature reaction modeling
- **Materials Science**: Phase equilibria and materials stability predictions
- **Environmental Science**: Modeling of atmospheric processes and pollution chemistry

## Comprehensive Logging

The project includes comprehensive logging features to assist with debugging and monitoring:

- **Log Directory**: All logs are stored in a configurable directory (default: `logs/`)
- **Log Levels**: Configurable log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- **Rotation**: Log files are automatically rotated when they reach a specified size
- **System Information**: Logs include detailed system and environment information
- **Error Tracing**: Comprehensive error handling with full stack traces
- **Performance Metrics**: Timing information for major operations

### Configuration

Logging can be configured in the YAML configuration files:

```yaml
# Logging Configuration in Species.yaml or EquilibriumCalculation.yaml
logging:
  directory: "logs"  # Directory for log files
  level: "INFO"      # Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
  file_prefix: "thermo_"  # Prefix for log files
  console_output: true  # Also output logs to console
  max_file_size_mb: 10  # Maximum log file size in MB
  backup_count: 5       # Number of backup log files to keep
  include_timestamps: true  # Include timestamps in log messages
```

## Reproducible Environments

To ensure reproducible results across different systems, we provide two methods:

### Requirements File

The `requirements.txt` file specifies exact versions of all dependencies:

```
cantera==2.6.0
numpy==1.24.3
pandas==2.0.3
matplotlib==3.7.2
PyYAML==6.0
requests==2.31.0
scipy==1.10.1
```

Install with:
```bash
pip install -r requirements.txt
```

### Docker Container

For maximum reproducibility, use our Docker container which provides an isolated environment with all dependencies pre-installed:

1. Build the Docker image:
```bash
docker build -t thermodynamics .
```

2. Run the complete workflow with a single command:
```bash
docker run -it --rm -v $(pwd):/app thermodynamics run-all
```

3. Run individual calculations:
```bash
# Generate thermodynamic data only
docker run -it --rm -v $(pwd):/app thermodynamics python thermo_generator.py

# Calculate equilibrium concentrations only
docker run -it --rm -v $(pwd):/app thermodynamics python EquilibriumCalculation.py
```

4. Interactive shell in the container:
```bash
docker run -it --rm -v $(pwd):/app thermodynamics bash
```

5. Access logs and results:
```bash
# View log files
docker run -it --rm -v $(pwd):/app thermodynamics ls -l logs/

# View results
docker run -it --rm -v $(pwd):/app thermodynamics ls -l results/
```

#### Docker Container Features:

- **Persistent Volumes**: The container uses persistent volume mounts for `/app/logs` and `/app/results` to preserve data between runs
- **Entrypoint Script**: Includes a convenience `run-all` command to execute the complete workflow
- **Debugging Tools**: Pre-installed with `vim` and `less` for log inspection
- **Proper Permissions**: Directories have write permissions to avoid issues with mounted volumes
- **UTF-8 Support**: Environment variables ensure proper character encoding for logs
- **Compatibility**: Works on Linux, macOS, and Windows (via Docker Desktop)

The container mounts your current directory (`$(pwd)`) to `/app` inside the container, so all outputs are saved to your local filesystem.

## Future Work

- Multi-pressure calculations and pressure-temperature phase diagrams
- Reaction pathway analysis for key chemical transformations
- Species sensitivity analysis to identify key components
- Integration with atmospheric transport models
- Web-based visualization interface for interactive exploration

## License

MIT License

Copyright (c) 2025 David Lary (davidlary@me.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.