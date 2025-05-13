FROM condaforge/mambaforge:latest

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update -y && \
    apt-get install -y \
    build-essential \    
    && apt-get install -y --no-install-recommends git \
    && apt-get purge -y --auto-remove \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get update -y \
    && apt-get install -y \
    wget \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    gnupg2 \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libgfortran5 \
    libcurl4-openssl-dev \
    libpq-dev \
    gfortran \
    libmagick++-dev \
    libssl-dev \
    libblas-dev \
    liblapack-dev \
    libnetcdf-dev \
    netcdf-bin \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    zlib1g-dev \
    libreadline-dev \
    libgsl-dev \
    libxt-dev

RUN LIB_PATH=$(ldconfig -p | grep libgfortran | head -n 1 | awk '{print $NF}') \
    && LIB_DIR=$(dirname $LIB_PATH) \
    && export LD_LIBRARY_PATH=$LIB_DIR \
    && ln -s $LIB_PATH $LIB_DIR/libgfortran.so

# Clear APT cache and lists to ensure a fresh pull
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add the R CRAN repository for Ubuntu
RUN apt-get update -y && apt-get install -y --no-install-recommends software-properties-common \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends r-base

RUN R -q -e "install.packages('BiocManager', dependencies=TRUE)" && \
    R -q -e "options(repos = c(CRAN_mirror = 'https://cloud.r-project.org', BioC_mirror = 'https://bioconductor.org'))" && \
    R -q -e "options(warn=2); BiocManager::install(version = '3.21')"

RUN R -q -e "options(warn=2); install.packages('doRNG')"
RUN R -q -e "options(warn=2); BiocManager::install('GENIE3', dependencies=TRUE, force=TRUE)"


# Set working directory
WORKDIR /app

# Copy environment file
COPY environment_prod.yml /app/environment.yml

# Create conda environment
RUN mamba env create -f /app/environment.yml && \
    mamba clean -afy

# Copy application files
COPY . /app/

# Run the application
CMD ["mamba", "run", "-n", "grn", "python", "-m", "app.main"]
