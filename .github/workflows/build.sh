#!/bin/bash

# TOGA2 Build Script for CI
# This script replicates the build process and can be tested locally

set -e  # Exit on any error

echo "=============================="
echo "Starting TOGA2 Build Process"
echo "=============================="

# Configuration
VENV_NAME="toga2"
PYTHON_VERSION="3.10"

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to print section headers
print_section() {
    echo "=============================="
    echo "$1"
    echo "=============================="
}

# 1. Install system dependencies
print_section "Installing System Dependencies"
if command_exists apt-get; then
    sudo apt-get update
    sudo apt-get install -y \
        gcc \
        g++ \
        make \
        wget \
        curl \
        git \
        bash \
        cargo \
        rustc \
        pkg-config \
        libhdf5-dev \
        python3-dev \
        python3-pip \
        python3-venv
elif command_exists yum; then
    sudo yum install -y \
        gcc \
        gcc-c++ \
        make \
        wget \
        curl \
        git \
        bash \
        cargo \
        rustc \
        pkgconfig \
        hdf5-devel \
        python3-devel \
        python3-pip
else
    echo "Warning: Could not detect package manager. Please ensure system dependencies are installed."
fi

# 2. Create and activate virtual environment
print_section "Setting up Python Virtual Environment"
if [ ! -d "$VENV_NAME" ]; then
    python3 -m venv "$VENV_NAME"
    echo "Created virtual environment: $VENV_NAME"
fi

# Activate virtual environment
source "$VENV_NAME/bin/activate"
echo "Activated virtual environment: $VENV_NAME"

# 3. Upgrade pip and install Python packages
print_section "Installing Python Packages"
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# 4. Make executable permissions
print_section "Setting Executable Permissions"
chmod +x check_dependencies.py
chmod +x toga2.py
chmod +x src/python/*.py
chmod +x src/python/modules/*.py

# 5. Check dependencies
print_section "Checking Dependencies"
python check_dependencies.py essentials
python check_dependencies.py managers

# 6. Build C components
print_section "Building C Components"
mkdir -p bin
if [ "$(uname -m)" = "aarch64" ]; then
    CFLAGS="-Wall -Wextra -O2 -g -std=c99 -arch arm64"
else
    CFLAGS="-Wall -Wextra -O2 -g -std=c99"
fi
gcc $CFLAGS -o bin/chain_filter_by_id src/c/chain_filter_by_id.c
gcc $CFLAGS -fPIC -shared -o src/python/modules/chain_coords_converter_slib.so src/c/chain_coords_converter_slib.c
gcc $CFLAGS -fPIC -shared -o src/python/modules/chain_bst_lib.so src/c/chain_bst_lib.c

# 7. Download UCSC binaries
print_section "Downloading UCSC Binaries"
mkdir -p bin
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/faToTwoBit
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/wigToBigWig
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigWigToWig
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigBedToBed
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedToBigBed
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/ixIxx
chmod +x bin/*

# 8. Build Cython components
print_section "Building Cython Components"
mkdir -p bin
python setup.py build_ext --build-lib=src

# 9. Build Rust components
print_section "Building Rust Components"
cd src/rust && cargo build --release
cd ../..
cd bed2gtf && cargo build --release
cd ..

# 10. Build CESAR2.0
print_section "Building CESAR2.0"
cd CESAR2.0 && make
cd ..

# 11. Install postoga rustools
print_section "Installing Postoga Rustools"
cd postoga/rustools
maturin develop --release
cd ../..

# 12. Train models
print_section "Training Models"
python src/python/train_model.py

# 13. Final integration test
print_section "Final Integration Test"
./toga2.py -h

# 14. Verify build artifacts
print_section "Verifying Build Artifacts"
echo "Checking if binaries exist:"
ls -la bin/
echo "Checking if shared libraries exist:"
ls -la src/python/modules/*.so
echo "Checking if Rust binaries exist:"
ls -la src/rust/target/release/
ls -la bed2gtf/target/release/

print_section "Build Completed Successfully!"
echo "TOGA2 has been built and tested successfully."
echo "Virtual environment '$VENV_NAME' is active and contains all dependencies."
