FROM python:3.10-slim

LABEL maintainer="David Lary <davidlary@me.com>"
LABEL description="Atmospheric Thermodynamic Equilibrium Calculator"
LABEL version="1.1.0"
LABEL git_repository="github.com/davidlary/Thermodynamics"

# Install system dependencies required for Cantera
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    swig \
    libboost-dev \
    # Tools for debugging and logs
    vim \
    less \
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

# Create directories for outputs and logs with proper permissions
RUN mkdir -p results/plots results/summary logs && \
    chmod -R 777 results logs

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=utf-8
ENV PYTHONDONTWRITEBYTECODE=1

# Add volume mounts for persistent data
VOLUME ["/app/logs", "/app/results"]

# Add entrypoint script for convenience
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

# Default command
CMD ["bash"]

# Usage instructions as comment
# To generate thermodynamic data:
#   docker run -it --rm -v $(pwd):/app thermodynamics python thermo_generator.py
# To calculate equilibrium concentrations:
#   docker run -it --rm -v $(pwd):/app thermodynamics python EquilibriumCalculation.py
# To run the complete workflow:
#   docker run -it --rm -v $(pwd):/app thermodynamics run-all