# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
- This repository contains thermodynamics-related data, particularly chemical species listings.
- The data is stored in YAML format.

## Commands
- Validate YAML files: `yamllint Species.yaml`
- Process YAML data: Use Python with PyYAML (`import yaml; with open('Species.yaml') as f: data = yaml.safe_load(f)`)
- Run thermodynamic calculations: No specific command available yet

## Code Style
- **YAML Formatting**: Use 2-space indentation for YAML files
- **Species Names**: Follow standard chemical nomenclature
- **List Structure**: Maintain hierarchical structure with species as list items
- **Comments**: Add descriptive comments above major sections 
- **Line Length**: Keep lines reasonably short (under 80 characters)
- **Validation**: Check chemical formulas for accuracy before adding
- **Grouping**: Consider grouping related species together with comments

## Adding New Species
- Verify chemical formula correctness before adding
- Add to appropriate section (neutral molecules, ions, etc.)
- Maintain alphabetical order within logical groupings