# Atmospheric Thermodynamic Equilibrium Calculator

This project calculates thermodynamic equilibrium concentrations of chemical species as a function of temperature and pressure using the Cantera library, with a special focus on accuracy in the cold temperature range (100-400K) for atmospheric research applications.

## Features

- Reads chemical species from a YAML file
- Generates NASA-9 polynomial thermodynamic data for Cantera
- Prioritized thermodynamic data lookup: Burcat → CEA → NASA → NIST
- Local caching of thermodynamic data for improved performance
- Enhanced temperature range support from 100K to 6000K
- Optimized accuracy in cold temperature range (100-400K) for atmospheric research
- Robust temperature range transition continuity verification and optimization
- Support for all relevant chemical elements and charged species
- No placeholder data - only real thermodynamic data from validated sources

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
3. Downloads and caches the complete Burcat thermodynamic database
4. For each species, retrieves thermodynamic data in this priority order:
   - Burcat's Thermodynamic Database (cached locally, updated weekly)
   - NASA CEA Database
   - NASA Database
   - NIST Database
5. Locally caches all retrieved data to speed up future runs
6. Generates a Cantera-compatible YAML file with the complete thermodynamic data
7. Includes detailed source attribution for each species in the output
8. Creates a data sources report summarizing which database was used for each species

### Thermodynamic Data Sources

The system prioritizes data sources in this order:

1. **Burcat Database**: The most comprehensive and up-to-date source of thermodynamic data, maintained by Alexander Burcat at the Technion - Israel Institute of Technology.
2. **NASA CEA Database**: Chemical Equilibrium with Applications database from NASA.
3. **NASA Thermodynamic Database**: General NASA thermodynamic data.
4. **NIST Chemistry WebBook**: National Institute of Standards and Technology data.

The entire Burcat database is downloaded and cached locally, with automatic version checking to ensure the data remains current without unnecessary downloads.

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

## Applications

This tool is particularly useful for:

- **Atmospheric Research**: Enhanced accuracy in the 100-400K range makes it ideal for stratospheric and mesospheric chemistry studies
- **Chemical Engineering**: Accurate prediction of equilibrium concentrations across wide temperature ranges
- **Combustion Analysis**: High-temperature reaction modeling
- **Materials Science**: Phase equilibria and materials stability predictions
- **Environmental Science**: Modeling of atmospheric processes and pollution chemistry

## Future Work

- Advanced equilibrium concentration calculations
- Temperature and pressure dependency visualization
- Reaction pathway analysis
- Graphical interface for atmospheric research applications
- Integration with atmospheric transport models

## License

MIT