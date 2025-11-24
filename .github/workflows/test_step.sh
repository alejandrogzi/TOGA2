#!/bin/bash

# Individual test scripts for TOGA2 CI based on Makefile targets

set -e

echo "=============================="
echo "Running $1"
echo "=============================="

case "$1" in
    "check_shell")
        if [ "$(echo $0)" != "-bash" ]; then
            echo "ERROR: TOGA2 currently only supports bash as operating shell"
            exit 1
        else
            echo "bash has been found to be a default current shell; shell type check successfully passed"
        fi
        ;;
    
    "check_essentials")
        python3 check_dependencies.py essentials
        ;;
    
    "check_managers")
        python3 check_dependencies.py managers
        ;;
    
    "check_python")
        python3 check_dependencies.py python --installation_mode
        ;;
    
    "check_third_party")
        python3 check_dependencies.py third_party
        ;;
    
    "chmod")
        chmod +x check_dependencies.py
        chmod +x toga2.py
        for exec in cesar_exec.py cesar_preprocess.py classify_chains.py feature_extractor.py fine_orthology_resolver.py get_contig_sizes.py predict_with_spliceai.py train_model.py; do
            chmod +x src/python/$exec
        done
        for exec in chain_bst_index.py get_names_from_bed.py; do
            chmod +x src/python/modules/$exec
        done
        ;;
    
    "install_binaries")
        mkdir -p bin
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/faToTwoBit
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/wigToBigWig
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigWigToWig
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigBedToBed
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedToBigBed
        wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/ixIxx
        chmod +x bin/*
        ;;
    
    "install_python")
        python3 -m pip install -r requirements.txt
        ;;
    
    "install_third_party")
        python3 check_dependencies.py install_third_party
        ;;
    
    "build_c")
        if [ "$(uname -m)" = "aarch64" ]; then
            CFLAGS="-Wall -Wextra -O2 -g -std=c99 -arch arm64"
        else
            CFLAGS="-Wall -Wextra -O2 -g -std=c99"
        fi
        mkdir -p bin && gcc $CFLAGS -o bin/chain_filter_by_id src/c/chain_filter_by_id.c
        gcc $CFLAGS -fPIC -shared -o src/python/modules/chain_coords_converter_slib.so src/c/chain_coords_converter_slib.c
        gcc $CFLAGS -fPIC -shared -o src/python/modules/chain_bst_lib.so src/c/chain_bst_lib.c
        ;;
    
    "build_cesar")
        cd CESAR2.0 && make
        cd ..
        ;;
    
    "build_cython")
        mkdir -p bin && python3 setup.py build_ext --build-lib=src
        ;;
    
    "build_rust")
        cd src/rust && cargo build --release
        cd ../..
        cd bed2gtf && cargo build --release
        cd ..
        ;;
    
    "install_postoga")
        source toga2/bin/activate
        cd postoga/rustools
        maturin develop --release
        cd ../..
        ;;
    
    "train_models")
        source toga2/bin/activate
        python3 src/python/train_model.py
        ;;
    
    "final_test")
        source toga2/bin/activate
        ./toga2.py -h
        ;;
    
    "verify_artifacts")
        echo "Checking if binaries exist:"
        ls -la bin/
        echo "Checking if shared libraries exist:"
        ls -la src/python/modules/*.so
        echo "Checking if Rust binaries exist:"
        ls -la src/rust/target/release/
        ls -la bed2gtf/target/release/
        ;;
    
    *)
        echo "Unknown target: $1"
        echo "Available targets: check_shell, check_essentials, check_managers, check_python, check_third_party, chmod, install_binaries, install_python, install_third_party, build_c, build_cesar, build_cython, build_rust, install_postoga, train_models, final_test, verify_artifacts"
        exit 1
        ;;
esac

echo "=============================="
echo "Completed $1 successfully"
echo "=============================="