# Thermodynamic Equilibrium Calculator

This project calculates thermodynamic equilibrium concentrations of chemical species as a function of temperature and pressure using the Cantera library.

## Features

- Reads chemical species from a YAML file
- Generates NASA-9 polynomial thermodynamic data for Cantera
- Prioritized thermodynamic data lookup: Burcat → CEA → NASA → NIST
- Local caching of thermodynamic data for improved performance
- Temperature range support from 100K to 6000K
- Support for all relevant chemical elements and charged species

## Requirements

- Python 3.6+
- Cantera
- PyYAML
- Requests

## Installation

```bash
pip install cantera pyyaml requests
```

## Usage

1. Define your species in `Species.yaml`:

```yaml
# List of chemical species
species:
  - N2
  - O2
  - H2O
  - CO2
  # Add more species as needed
```

2. Generate thermodynamic data:

```bash
python thermo_generator.py
```

This will create a `Thermodynamics.yaml` file with the NASA-9 polynomial data ready for use in Cantera.

3. Use the data for equilibrium calculations (future implementation).

## How It Works

The system follows these steps:

1. Reads the list of species from `Species.yaml`
2. Extracts all unique elements from the chemical formulas
3. For each species, retrieves thermodynamic data in this priority order:
   - Burcat's Thermodynamic Database
   - NASA CEA Database
   - NASA Database
   - NIST Database
4. Locally caches all retrieved data to speed up future runs
5. Generates a Cantera-compatible YAML file with the complete thermodynamic data

## Future Work

- Equilibrium concentration calculations
- Temperature and pressure dependency analysis
- Reaction pathway visualization
- Web interface for calculations

## License

MIT