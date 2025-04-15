#!/usr/bin/env python3
"""
Equilibrium Concentration Calculator

This script calculates chemical equilibrium concentrations over a specified
temperature range using Cantera and the NASA-9 polynomial thermodynamic data.
Configuration parameters are read from EquilibriumCalculation.yaml.
"""

import os
import sys
import time
import yaml
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import cantera as ct
import traceback
from logging.handlers import RotatingFileHandler
import datetime
import platform

# Logger setup will be done after loading configuration
logger = logging.getLogger('equilibrium_calc')

def setup_logging(config: Dict) -> None:
    """
    Set up comprehensive logging based on configuration.
    
    Args:
        config: Configuration dictionary with logging settings
    """
    # Get logging configuration
    log_config = config.get('logging', {})
    log_dir = log_config.get('directory', 'logs')
    log_level_str = log_config.get('level', 'INFO').upper()
    file_prefix = log_config.get('file_prefix', 'equilibrium_')
    console_output = log_config.get('console_output', True)
    max_file_size = log_config.get('max_file_size_mb', 10) * 1024 * 1024  # Convert MB to bytes
    backup_count = log_config.get('backup_count', 5)
    include_timestamps = log_config.get('include_timestamps', True)
    
    # Map string log level to logging constant
    log_level_map = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR,
        'CRITICAL': logging.CRITICAL
    }
    log_level = log_level_map.get(log_level_str, logging.INFO)
    
    # Create log directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)
    
    # Define log format
    if include_timestamps:
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    else:
        log_format = '%(name)s - %(levelname)s - %(message)s'
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(log_format)
    
    # Set up file handler with rotation
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f"{file_prefix}{timestamp}.log")
    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=max_file_size,
        backupCount=backup_count
    )
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
    
    # Set up console handler if enabled
    if console_output:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
    
    # Log system information
    logger.info("=" * 80)
    logger.info("Equilibrium Concentration Calculator - Started")
    logger.info("=" * 80)
    logger.info(f"Date and Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Python Version: {sys.version}")
    logger.info(f"Platform: {platform.platform()}")
    logger.info(f"Cantera Version: {ct.__version__}")
    logger.info(f"Numpy Version: {np.__version__}")
    logger.info(f"Pandas Version: {pd.__version__}")
    logger.info(f"Matplotlib Version: {plt.matplotlib.__version__}")
    logger.info(f"Log Level: {log_level_str}")
    logger.info(f"Log File: {log_file}")
    logger.info("=" * 80)
    
    return

def load_config(config_file: str = 'EquilibriumCalculation.yaml') -> Dict:
    """
    Load configuration from YAML file.
    
    Args:
        config_file: Path to YAML configuration file
        
    Returns:
        Dictionary containing configuration parameters
    """
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded configuration from {config_file}")
        return config
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        sys.exit(1)

def validate_config(config: Dict) -> Dict:
    """
    Validate configuration parameters and set defaults.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Validated configuration dictionary
    """
    # Required sections
    required_sections = ['input', 'calculation', 'output', 'initial_mixture']
    for section in required_sections:
        if section not in config:
            logger.error(f"Missing required section: {section}")
            sys.exit(1)
    
    # Input section
    if 'thermo_data_file' not in config['input']:
        config['input']['thermo_data_file'] = "Thermodynamics.yaml"
    
    # Add use_cantera_database option if not present - default to FALSE to use our custom data
    if 'use_cantera_database' not in config['input']:
        config['input']['use_cantera_database'] = False
        logger.info("Using custom thermodynamic data from Thermodynamics.yaml by default")
        
    # Check if thermodynamic data file exists
    thermo_file = config['input']['thermo_data_file']
    if not os.path.exists(thermo_file):
        logger.error(f"Thermodynamic data file not found: {thermo_file}")
        if not config['input']['use_cantera_database']:
            logger.warning("Critical error: Custom thermodynamic data file not found")
            logger.warning("Please run thermo_generator.py first to create Thermodynamics.yaml")
            sys.exit(1)
        else:
            logger.warning("Using Cantera's built-in database since specified thermodynamic file not found")
    
    # Calculation section
    calc = config['calculation']
    if 'temperature_range' not in calc:
        calc['temperature_range'] = [100, 6000]
    if 'temperature_points' not in calc:
        calc['temperature_points'] = 500
    if 'log_temperature' not in calc:
        calc['log_temperature'] = True
    if 'pressure' not in calc:
        calc['pressure'] = 101325  # 1 atm in Pa
    if 'solver_tolerance' not in calc:
        calc['solver_tolerance'] = 1.0e-8
    if 'max_iterations' not in calc:
        calc['max_iterations'] = 500
    
    # Output section
    out = config['output']
    if 'directory' not in out:
        out['directory'] = "results"
    if 'csv_file' not in out:
        out['csv_file'] = "equilibrium_results.csv"
    if 'save_csv' not in out:
        out['save_csv'] = True
    if 'plot_results' not in out:
        out['plot_results'] = out.get('create_plots', True)
    if 'create_plots' not in out:
        out['create_plots'] = out.get('plot_results', True)
    if 'plot_filename' not in out:
        out['plot_filename'] = "equilibrium_plot.png"
    if 'plot_format' not in out:
        out['plot_format'] = "png"
    if 'plots_per_page' not in out:
        out['plots_per_page'] = 24
    if 'plots_per_row' not in out:
        out['plots_per_row'] = 4
    if 'log_scale' not in out:
        out['log_scale'] = True
    if 'plot_scale' not in out:
        out['plot_scale'] = "log"
    if 'plot_type' not in out:
        out['plot_type'] = "concentration"
    if 'focus_species' not in out:
        out['focus_species'] = out.get('species_to_plot', [])
    if 'species_to_plot' not in out:
        out['species_to_plot'] = out.get('focus_species', [])
    if 'exclude_species' not in out:
        out['exclude_species'] = []
    if 'concentration_threshold' not in out:
        out['concentration_threshold'] = 1.0e-12
    if 'log_max_concentrations' not in out:
        out['log_max_concentrations'] = True
    if 'log_initial_mixture' not in out:
        out['log_initial_mixture'] = True
    
    # Logging section (add if not present)
    if 'logging' not in config:
        config['logging'] = {}
    
    log_config = config['logging']
    if 'directory' not in log_config:
        log_config['directory'] = "logs"
    if 'level' not in log_config:
        log_config['level'] = "INFO"
    if 'file_prefix' not in log_config:
        log_config['file_prefix'] = "equilibrium_"
    if 'console_output' not in log_config:
        log_config['console_output'] = True
    if 'max_file_size_mb' not in log_config:
        log_config['max_file_size_mb'] = 10
    if 'backup_count' not in log_config:
        log_config['backup_count'] = 5
    if 'include_timestamps' not in log_config:
        log_config['include_timestamps'] = True
    
    # Validate log level
    valid_log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    log_config['level'] = log_config['level'].upper()
    if log_config['level'] not in valid_log_levels:
        logger.warning(f"Invalid log level: {log_config['level']}. Using INFO instead.")
        log_config['level'] = "INFO"
    
    # Validate initial_mixture
    mixture = config['initial_mixture']
    total = sum(mixture.values())
    if not 0.99 <= total <= 1.01:
        logger.warning(f"Initial mixture mole fractions sum to {total}, not 1.0")
        # Normalize
        for species in mixture:
            mixture[species] /= total
        logger.info("Normalized initial mixture mole fractions to sum to 1.0")
    
    return config

def setup_output_directory(config: Dict) -> str:
    """
    Set up output directory for results.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Path to output directory
    """
    output_dir = config['output']['directory']
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories for plots if needed
    if config['output']['create_plots']:
        os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'summary'), exist_ok=True)
    
    logger.info(f"Output will be saved to {output_dir}")
    return output_dir

def create_temperature_array(config: Dict) -> np.ndarray:
    """
    Create temperature array for calculations.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        NumPy array of temperatures
    """
    t_min, t_max = config['calculation']['temperature_range']
    n_points = config['calculation']['temperature_points']
    log_spacing = config['calculation'].get('log_temperature', True)
    
    if log_spacing:
        # Use logarithmic spacing for better resolution at low temperatures
        # which is important for atmospheric applications
        t_array = np.logspace(np.log10(t_min), np.log10(t_max), n_points)
        spacing_type = "logarithmic"
    else:
        # Use linear spacing
        t_array = np.linspace(t_min, t_max, n_points)
        spacing_type = "linear"
    
    logger.info(f"Created temperature array with {n_points} points ({spacing_type} spacing) from {t_min} K to {t_max} K")
    return t_array

def create_simple_gas_model(config: Dict) -> str:
    """
    Create a simple gas model file using known species from Cantera's built-in database.
    
    Instead of trying to convert our complex YAML file, we'll use Cantera's
    built-in species database and just create a list of the species we need.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Path to created XML file
    """
    # Get list of species from initial mixture
    species_list = list(config['initial_mixture'].keys())
    
    # Load all species from Species.yaml if available
    try:
        with open('Species.yaml', 'r') as f:
            species_data = yaml.safe_load(f)
            if 'species' in species_data and isinstance(species_data['species'], list):
                # Add all species from the file to ensure we have a complete set
                species_list.extend([s for s in species_data['species'] if s not in species_list])
                logger.info(f"Added {len(species_data['species'])} species from Species.yaml")
    except Exception as e:
        logger.warning(f"Could not load Species.yaml: {e}")
    
    # Create a temporary XML file
    xml_file = "simple_gas.xml"
    
    # Check for charged species
    has_charged_species = any('-' in s or '+' in s or s == 'e-' for s in species_list)
    
    # Create a simple XML structure
    xml_content = """<?xml version="1.0"?>
<ctml>
  <validate reactions="yes" species="yes"/>

  <phase dim="3" id="gas">
    <elementArray datasrc="elements.xml">E H C N O Ar He Ne Xe S</elementArray>
    <speciesArray datasrc="#species_data">"""
    
    # Add species
    for species in species_list:
        xml_content += f" {species}"
    
    xml_content += """</speciesArray>
    <reactionArray datasrc="#reaction_data"/>
    <state>
      <temperature units="K">300.0</temperature>
      <pressure units="Pa">101325.0</pressure>
    </state>
    <thermo model="IdealGas"/>
    <kinetics model="GasKinetics"/>
    <transport model="None"/>"""
    
    # Add electron energy mode if we have charged species
    if has_charged_species:
        xml_content += """
    <electron_energy enabled="yes"/>"""
    
    xml_content += """
  </phase>

  <speciesData id="species_data">
    <!-- Simplified species data - will use Cantera's built-in data for known species -->
  </speciesData>

  <reactionData id="reaction_data"/>
</ctml>
"""
    
    # Write XML file
    with open(xml_file, 'w') as f:
        f.write(xml_content)
    
    logger.info(f"Created simplified gas model in {xml_file} with {len(species_list)} species")
    
    # Print warning about potential missing species
    logger.warning("Using Cantera's built-in database. Some species may not be available.")
    logger.info("If equilibrium calculation fails, check if all required species are in Cantera's database.")
    
    return xml_file

def convert_yaml_to_cantera_yaml(config: Dict) -> str:
    """
    Convert our custom YAML thermodynamic data to Cantera-compatible YAML format.
    
    This function converts our NASA-9 format thermodynamic data in YAML to 
    a format that Cantera can read (Cantera-compatible YAML).
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Path to created Cantera-compatible YAML file
    """
    thermo_file = config['input']['thermo_data_file']
    logger.info(f"Converting thermodynamic data from {thermo_file} to Cantera format")
    
    try:
        # Load our custom thermodynamic data
        with open(thermo_file, 'r') as f:
            thermo_data = yaml.safe_load(f)
        
        # Create Cantera YAML file path
        cantera_yaml_file = "converted_thermo.yaml"
        
        # Create Cantera-compatible YAML structure
        cantera_yaml = {
            "description": "Converted from custom thermodynamic data",
            "generator": "EquilibriumCalculation.py",
            "cantera-version": "3.0",
            "date": datetime.datetime.now().strftime("%Y-%m-%d"),
            "units": {"length": "cm", "time": "s", "quantity": "mol", "energy": "cal"},
            "phases": [
                {
                    "name": "gas",
                    "thermo": "ideal-gas",
                    "elements": ["H", "C", "N", "O", "Ar", "He", "Ne", "Xe", "S", "E"],
                    "species": [],
                    "state": {"T": 300.0, "P": 101325.0}
                }
            ],
            "species": []
        }
        
        # Process each species
        if 'species' in thermo_data:
            species_list = []
            for species in thermo_data['species']:
                name = species['name']
                logger.info(f"Processing species: {name}")
                
                # Add to phase species list
                cantera_yaml["phases"][0]["species"].append(name)
                species_list.append(name)
                
                # Create species entry
                species_entry = {
                    "name": name,
                    "composition": species.get('composition', {})
                }
                
                # Process NASA polynomial data
                if 'thermo' in species and 'temperature-ranges' in species['thermo']:
                    # Convert to Cantera's NASA polynomial format
                    nasa_polys = []
                    
                    for temp_range in species['thermo']['temperature-ranges']:
                        t_min = temp_range['T-min']
                        t_max = temp_range['T-max']
                        coeffs = temp_range['coefficients']
                        
                        # Ensure we have 7 coefficients (Cantera's NASA format is 7 coefficients)
                        nasa7_coeffs = coeffs[:7] if len(coeffs) >= 7 else coeffs + [0.0] * (7 - len(coeffs))
                        
                        nasa_poly = {
                            "T_range": [t_min, t_max],
                            "coeffs": nasa7_coeffs
                        }
                        nasa_polys.append(nasa_poly)
                    
                    # Sort polynomials by temperature range
                    nasa_polys = sorted(nasa_polys, key=lambda x: x["T_range"][0])
                    
                    # Add thermo data to species entry
                    species_entry["thermo"] = {
                        "model": "NASA7",
                        "reference-pressure": 101325.0,  # 1 atm in Pa
                        "temperature-ranges": [poly["T_range"] for poly in nasa_polys],
                        "data": [poly["coeffs"] for poly in nasa_polys]
                    }
                
                # Add to species list
                cantera_yaml["species"].append(species_entry)
            
            logger.info(f"Processed {len(species_list)} species")
        
        # Write Cantera YAML file
        with open(cantera_yaml_file, 'w') as f:
            yaml.dump(cantera_yaml, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Successfully converted thermodynamic data to {cantera_yaml_file}")
        return cantera_yaml_file
        
    except Exception as e:
        logger.error(f"Error converting to Cantera YAML: {e}")
        logger.error(traceback.format_exc())
        logger.error("Falling back to simplified gas model with Cantera's built-in database")
        return create_simple_gas_model(config)

def setup_cantera_gas(config: Dict) -> ct.Solution:
    """
    Set up Cantera gas mixture with thermodynamic data.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Cantera Solution object
    """
    try:
        # Check if we should use our custom thermodynamic data or Cantera's built-in database
        use_cantera_database = config.get('input', {}).get('use_cantera_database', False)
        
        if not use_cantera_database:
            # Convert our custom thermodynamic data to a format Cantera can read
            logger.info("Using custom thermodynamic data from Thermodynamics.yaml")
            
            try:
                # Try direct loading first (in case format is already compatible)
                thermo_file = config['input']['thermo_data_file']
                try:
                    gas = ct.Solution(thermo_file)
                    logger.info(f"Successfully loaded thermodynamic data directly from {thermo_file}")
                except Exception as direct_error:
                    logger.info(f"Could not load thermodynamic data directly: {direct_error}")
                    logger.info("Converting YAML format to Cantera-compatible format...")
                    
                    # Convert our YAML to Cantera-compatible YAML format
                    cantera_yaml_file = convert_yaml_to_cantera_yaml(config)
                    gas = ct.Solution(cantera_yaml_file)
                    logger.info(f"Successfully loaded thermodynamic data from converted file {cantera_yaml_file}")
            except Exception as e:
                logger.error(f"Error loading custom thermodynamic data: {e}")
                logger.warning("IMPORTANT: Falling back to Cantera's built-in database as a last resort")
                xml_file = create_simple_gas_model(config)
                gas = ct.Solution(xml_file)
                logger.warning("Using Cantera's built-in database instead of our custom thermodynamic data")
        else:
            # Explicitly using Cantera's built-in database (this should not happen based on config)
            logger.warning("IMPORTANT NOTE: Using Cantera's built-in database as specified in configuration")
            xml_file = create_simple_gas_model(config)
            gas = ct.Solution(xml_file)
        
        logger.info(f"Created Cantera gas mixture with {gas.n_species} species")
        
        # Set initial state
        P = config['calculation']['pressure']
        T = 298.15  # Standard temperature for initial state
        
        # Convert mole fractions to string format for Cantera
        X = config['initial_mixture']
        # Filter to only include species in the gas model
        X_filtered = {k: v for k, v in X.items() if k in gas.species_names}
        
        # Log information about missing species
        missing_species = [k for k in X.keys() if k not in gas.species_names]
        if missing_species:
            logger.warning(f"The following species from the initial mixture were not found in the gas model:")
            for species in missing_species:
                logger.warning(f"  - {species}")
        
        # Normalize the filtered mole fractions
        total = sum(X_filtered.values())
        if total > 0:
            X_filtered = {k: v/total for k, v in X_filtered.items()}
        
        # Convert to string format
        X_str = ", ".join([f"{k}:{v}" for k, v in X_filtered.items()])
        
        # Set state
        gas.TPX = T, P, X_str
        
        # Log initial mixture if requested
        if config['output']['log_initial_mixture']:
            logger.info(f"Initial mixture composition (mole fractions):")
            for species, mole_frac in sorted(X.items()):
                if species in gas.species_names:
                    logger.info(f"  {species}: {mole_frac:.6e}")
                else:
                    logger.warning(f"  {species}: not found in the gas model")
        
        # Check if we have a valid mixture
        if not X_filtered:
            logger.error("None of the specified species were found in the gas model")
            logger.info("Available species: " + ", ".join(gas.species_names))
            sys.exit(1)
            
        return gas
    
    except Exception as e:
        logger.error(f"Error setting up Cantera gas mixture: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)

def calculate_equilibrium(
    gas: ct.Solution, 
    temperatures: np.ndarray, 
    pressure: float,
    solver_tolerance: float = 1.0e-8,
    max_iterations: int = 500
) -> pd.DataFrame:
    """
    Calculate equilibrium concentrations over temperature range.
    
    Args:
        gas: Cantera Solution object
        temperatures: Array of temperatures (K)
        pressure: Pressure (Pa)
        solver_tolerance: Solver tolerance
        max_iterations: Maximum iterations for equilibrium solver
        
    Returns:
        DataFrame with equilibrium concentrations
    """
    # Store results
    n_temps = len(temperatures)
    n_species = gas.n_species
    species_names = gas.species_names
    
    # Create empty DataFrame
    columns = ['Temperature', 'Pressure'] + species_names
    result_df = pd.DataFrame(columns=columns)
    
    # Calculate equilibrium at each temperature
    for i, T in enumerate(temperatures):
        try:
            # Set initial state
            gas.TP = T, pressure
            
            # Calculate equilibrium
            gas.equilibrate('TP', solver='gibbs', max_steps=max_iterations, rtol=solver_tolerance)
            
            # Get mole fractions
            X = gas.X
            
            # Store results
            row = {'Temperature': T, 'Pressure': pressure}
            for j, species in enumerate(species_names):
                row[species] = X[j]
                
            result_df = pd.concat([result_df, pd.DataFrame([row])], ignore_index=True)
            
            # Log progress
            if (i + 1) % 50 == 0 or i == 0 or i == n_temps - 1:
                logger.info(f"Calculated equilibrium at T = {T:.2f} K ({i+1}/{n_temps})")
                
        except Exception as e:
            logger.error(f"Error calculating equilibrium at T = {T:.2f} K: {e}")
            # Add row with NaNs for this temperature
            row = {'Temperature': T, 'Pressure': pressure}
            for species in species_names:
                row[species] = np.nan
            result_df = pd.concat([result_df, pd.DataFrame([row])], ignore_index=True)
    
    return result_df

def save_results_to_csv(result_df: pd.DataFrame, config: Dict, output_dir: str) -> str:
    """
    Save results to CSV file.
    
    Args:
        result_df: DataFrame with equilibrium results
        config: Configuration dictionary
        output_dir: Output directory
        
    Returns:
        Path to saved CSV file
    """
    csv_file = os.path.join(output_dir, config['output']['csv_file'])
    result_df.to_csv(csv_file, index=False)
    logger.info(f"Saved results to {csv_file}")
    return csv_file

def create_species_plots(
    result_df: pd.DataFrame, 
    config: Dict, 
    output_dir: str
) -> List[str]:
    """
    Create concentration vs temperature plots for each species.
    
    Args:
        result_df: DataFrame with equilibrium results
        config: Configuration dictionary
        output_dir: Output directory
        
    Returns:
        List of plot file paths
    """
    if not config['output']['create_plots']:
        return []
    
    plots_dir = os.path.join(output_dir, 'plots')
    plot_format = config['output']['plot_format']
    log_scale = config['output']['log_scale']
    threshold = config['output']['concentration_threshold']
    focus_species = config['output']['focus_species']
    exclude_species = config['output']['exclude_species']
    
    # Get temperature and species columns
    temperature = result_df['Temperature']
    species_columns = [col for col in result_df.columns 
                    if col not in ['Temperature', 'Pressure']]
    
    # Filter species based on configuration
    if focus_species:
        # Only include specified species
        species_columns = [s for s in species_columns if s in focus_species]
    
    # Exclude specified species
    species_columns = [s for s in species_columns if s not in exclude_species]
    
    # Filter species based on maximum concentration
    if threshold > 0:
        species_columns = [s for s in species_columns 
                        if result_df[s].max() >= threshold]
    
    plot_files = []
    
    # Create individual plots for each species
    for species in species_columns:
        try:
            plt.figure(figsize=(10, 6))
            
            # Plot concentration vs temperature
            plt.plot(temperature, result_df[species], 'b-', linewidth=2)
            
            # Set axes and labels
            plt.xlabel('Temperature (K)')
            plt.ylabel('Mole Fraction')
            plt.title(f'{species} Equilibrium Concentration')
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            
            # Use log scale if specified
            if log_scale:
                plt.yscale('log')
                # Set reasonable y-limits
                y_min = max(result_df[species].min(), 1e-30)
                y_max = max(result_df[species].max() * 2, y_min * 10)
                plt.ylim(y_min, y_max)
            
            # Save plot
            plot_file = os.path.join(plots_dir, f'{species}.{plot_format}')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_files.append(plot_file)
            
        except Exception as e:
            logger.error(f"Error creating plot for {species}: {e}")
    
    logger.info(f"Created {len(plot_files)} species plots in {plots_dir}")
    return plot_files

def create_summary_plots(
    result_df: pd.DataFrame, 
    config: Dict, 
    output_dir: str
) -> List[str]:
    """
    Create summary plots with multiple species per page.
    
    Args:
        result_df: DataFrame with equilibrium results
        config: Configuration dictionary
        output_dir: Output directory
        
    Returns:
        List of plot file paths
    """
    if not config['output']['create_plots']:
        return []
    
    summary_dir = os.path.join(output_dir, 'summary')
    plot_format = config['output']['plot_format']
    log_scale = config['output']['log_scale']
    threshold = config['output']['concentration_threshold']
    plots_per_page = config['output']['plots_per_page']
    plots_per_row = config['output']['plots_per_row']
    focus_species = config['output']['focus_species']
    exclude_species = config['output']['exclude_species']
    
    # Get temperature and species columns
    temperature = result_df['Temperature']
    species_columns = [col for col in result_df.columns 
                    if col not in ['Temperature', 'Pressure']]
    
    # Filter species based on configuration
    if focus_species:
        # Only include specified species
        species_columns = [s for s in species_columns if s in focus_species]
    
    # Exclude specified species
    species_columns = [s for s in species_columns if s not in exclude_species]
    
    # Filter species based on maximum concentration
    if threshold > 0:
        species_columns = [s for s in species_columns 
                        if result_df[s].max() >= threshold]
    
    # Sort species by maximum concentration
    species_columns = sorted(species_columns, 
                          key=lambda s: result_df[s].max(), 
                          reverse=True)
    
    plot_files = []
    
    # Calculate number of pages needed
    n_species = len(species_columns)
    n_pages = (n_species + plots_per_page - 1) // plots_per_page
    
    for page in range(n_pages):
        # Get species for this page
        start_idx = page * plots_per_page
        end_idx = min((page + 1) * plots_per_page, n_species)
        page_species = species_columns[start_idx:end_idx]
        
        # Calculate grid dimensions
        n_rows = (len(page_species) + plots_per_row - 1) // plots_per_row
        
        # Create figure
        fig = plt.figure(figsize=(15, 3 * n_rows))
        gs = gridspec.GridSpec(n_rows, plots_per_row)
        
        # Create subplot for each species
        for i, species in enumerate(page_species):
            row = i // plots_per_row
            col = i % plots_per_row
            
            ax = fig.add_subplot(gs[row, col])
            
            # Plot concentration vs temperature
            ax.plot(temperature, result_df[species], 'b-', linewidth=1.5)
            
            # Set axes and labels
            ax.set_xlabel('Temperature (K)')
            ax.set_ylabel('Mole Fraction')
            ax.set_title(species)
            ax.grid(True, which='both', linestyle='--', alpha=0.5)
            
            # Use log scale if specified
            if log_scale:
                ax.set_yscale('log')
                # Set reasonable y-limits
                y_min = max(result_df[species].min(), 1e-30)
                y_max = max(result_df[species].max() * 2, y_min * 10)
                ax.set_ylim(y_min, y_max)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(summary_dir, f'summary_page{page+1}.{plot_format}')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        plot_files.append(plot_file)
    
    logger.info(f"Created {len(plot_files)} summary plots in {summary_dir}")
    return plot_files

def log_max_concentrations(result_df: pd.DataFrame, config: Dict) -> None:
    """
    Log maximum concentrations for each species.
    
    Args:
        result_df: DataFrame with equilibrium results
        config: Configuration dictionary
    """
    if not config['output']['log_max_concentrations']:
        return
    
    # Get species columns
    species_columns = [col for col in result_df.columns 
                    if col not in ['Temperature', 'Pressure']]
    
    # Calculate maximum concentration and corresponding temperature for each species
    max_conc = {}
    for species in species_columns:
        idx = result_df[species].idxmax()
        max_value = result_df[species].iloc[idx]
        temp = result_df['Temperature'].iloc[idx]
        max_conc[species] = (max_value, temp)
    
    # Sort by maximum concentration
    sorted_species = sorted(max_conc.items(), key=lambda x: x[1][0], reverse=True)
    
    # Log results
    logger.info("Maximum concentrations (mole fractions) for each species:")
    for species, (max_value, temp) in sorted_species:
        logger.info(f"  {species}: {max_value:.6e} at {temp:.2f} K")

def main():
    """Main function to run equilibrium calculations."""
    start_time = time.time()
    
    try:
        # Load and validate configuration
        logger.debug("Loading configuration...")
        config = load_config()
        
        # Set up logging based on configuration
        setup_logging(config)
        
        logger.debug("Validating configuration...")
        config = validate_config(config)
        logger.info("Configuration loaded and validated")
        
        # Log the configuration details at debug level
        logger.debug("Configuration details:")
        for section, settings in config.items():
            logger.debug(f"  {section}:")
            if isinstance(settings, dict):
                for key, value in settings.items():
                    logger.debug(f"    {key}: {value}")
            else:
                logger.debug(f"    {settings}")
        
        # Set up output directory
        logger.debug("Setting up output directory...")
        output_dir = setup_output_directory(config)
        
        # Create temperature array
        logger.debug("Creating temperature array...")
        temperatures = create_temperature_array(config)
        
        # Set up Cantera gas mixture
        logger.debug("Setting up Cantera gas mixture...")
        gas = setup_cantera_gas(config)
        
        # Calculate equilibrium concentrations
        pressure = config['calculation']['pressure']
        solver_tolerance = config['calculation']['solver_tolerance']
        max_iterations = config['calculation']['max_iterations']
        
        logger.info("Starting equilibrium calculations...")
        result_df = calculate_equilibrium(
            gas, 
            temperatures, 
            pressure,
            solver_tolerance,
            max_iterations
        )
        
        # Save results to CSV
        logger.debug("Saving results to CSV...")
        csv_file = save_results_to_csv(result_df, config, output_dir)
        
        # Create plots
        if config['output']['create_plots']:
            logger.debug("Creating species plots...")
            species_plots = create_species_plots(result_df, config, output_dir)
            
            logger.debug("Creating summary plots...")
            summary_plots = create_summary_plots(result_df, config, output_dir)
            
            logger.info(f"Created {len(species_plots)} species plots and {len(summary_plots)} summary plots")
        
        # Log maximum concentrations
        logger.debug("Logging maximum concentrations...")
        log_max_concentrations(result_df, config)
        
        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f"Equilibrium calculations completed in {elapsed_time:.2f} seconds")
        logger.info(f"Results saved to {output_dir}")
        logger.info("=" * 80)
        logger.info("Equilibrium Concentration Calculator - Completed Successfully")
        logger.info("=" * 80)
        
    except Exception as e:
        elapsed_time = time.time() - start_time
        logger.error("=" * 80)
        logger.error("ERROR: Equilibrium calculation failed!")
        logger.error(f"Elapsed time before error: {elapsed_time:.2f} seconds")
        logger.error(f"Error details: {str(e)}")
        logger.error("Traceback:")
        logger.error(traceback.format_exc())
        logger.error("=" * 80)
        sys.exit(1)

if __name__ == "__main__":
    # Set up basic logging until config is loaded
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    try:
        main()
    except Exception as e:
        logger.critical(f"Unhandled exception: {str(e)}")
        logger.critical(traceback.format_exc())
        sys.exit(1)