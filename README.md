# TOGA2
TOGA2: A faster, more versatile successor of Tool to infer Orthologs from Genome Alignments

> [!WARNING]  
TOGA2 is currently in early access phase. This means, certain TOGA2 features (`cookbook` and `from-congig` TOGA2 modes, automatic input preparation, etc.) and most of the documentation at TOGA2 wiki are currently missing and will be added in the following days. If you want to use the early access version and experience problems with installing or running code in this repository, please open an issue here or contact the TOGA2 team (yury.malovichko@senckenberg.de).

## Installation

### via `git clone`
> [!WARNING]  
> Local installation is under development and has not been thoroughly tested. As of now, we highly recommend running TOGA2 with Apptainer (see below).
> 
TOGA2 local installation (currently) relies on Make. The provided Makefile contains the following directives:
* `build` (default directive) checks for necessary third-party software and Python packages, compiles CESAR2 and TOGA2 modules;
* `install` fetches missing third-party software and installs missing/outdated Python packages 
```bash
git clone https://github.com/hillerlab/TOGA2
cd TOGA2
make
make install ## optional
```
#### Essential third-party software
`build` directive first checks for availability of the following commands in $PATH:
* awk
* cargo
* gcc
These utilities are necessary for further TOGA2 installation. If any of those are missing, the installation halts, prompting the user to install them manually

### via Apptainer
Container image definition file for Apptainer is provided with TOGA2 under `supply/apptainer.def`. To build the container, make sure you have Apptainer installed, then run the following commands:
```bash
git clone https://github.com/hillerlab/TOGA2
TMPDIR=${tmp_dir} apptainer build ${container.sif} supply/apptainer.def
```
where:
* `${tmp_dir}` is the path to the directory to store Nextflow temporary files in;
* `${container.sif}` is the path to save the resulting image to.

The resulting container has `toga2.py` as an entry point. If you run the following or similar command:
```bash
TMPDIR=${tmp_dir} apptainer run --bind ${bound_dir1},${bound_dir2} ${container.sif}
```
you should see the TOGA2 start menu.

## Test run
> [!WARNING]
> Currently not implemented


To ensure the proper TOGA2 installation/container build, run the following command:
```bash
toga2.py test
```
or its Apptainer alternative, respectively.
