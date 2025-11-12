# TOGA2
TOGA2: A faster, more versatile successor of Tool to infer Orthologs from Genome Alignments

> [!IMPORTANT]  
TOGA2 is currently in early access phase. This means, certain TOGA2 features and most of the documentation at TOGA2 wiki are currently missing and will be added in the following days. If you want to use the early access version and experience problems with installing or running code in this repository, please open an issue here or contact the TOGA2 team (yury.malovichko@senckenberg.de).

## Installation

### via `git clone`
> [!IMPORTANT]  
> Local installation is currently under development. If you encounter any problems with installation, please open a new Issue or contact TOGA2 developers
> 

```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
cd TOGA2
make
```
Since the Make directive also installs Python packages globally, you might want to create a dedicated virtual environment beforehand (Python version 3.9 or higher:
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
cd TOGA2
python3 -m venv toga2
source toga2/bin/activate
make
```

### via Apptainer
Container image definition file for Apptainer is provided with TOGA2 under `supply/apptainer.def`. To build the container, make sure you have Apptainer installed, then run the following commands:
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
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

>[!NOTE]
> The image provided in `supply/` directory contains the latest TOGA2 release, third-party software used for input preparation and TOGA2 annotation, and Nextflow for parallel process management. The container does <ins>not</ins> contain any Nextflow-compatible parallel job executor, and running parallel jobs from the container seems highly unlikely due to the code organization.  

## Running TOGA2
If activated without additional arguments (`toga2.py`), the following start screen is displayed in the user's terminal:
```
Usage: toga2.py [OPTIONS] COMMAND [ARGS]...

  MMP""MM""YMM   .g8""8q.     .g8"""bgd      db          `7MMF'`7MMF'
  P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
       MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
       MM     MM        MM MM             ,M  `MM          MM    MM  
       MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
       MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
     .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

  TOGA2 - Tool for Ortholog Inference from Genome Alignment

Options:
  -V, --version      Show the version and exit.
  -help, -h, --help  Show this message and exit.

Commands:
  cookbook       List example commands for 'run' mode
  from-config    Run TOGA2 pipeline with a configuration file
  postoga        Run postprocessing analysis with Postoga
  prepare-input  Prepare reference annotation files for TOGA2 input
  run            Run TOGA2 pipeline with command line arguments
  spliceai       Generate SpliceAI predictions for query assembly
  test           Test TOGA2 with companion dataset
```
Except for `toga2.py test`, invoking any of the listed commands without arguments also displays the help message. You can also invoke help message for TOGA2 or any of its commands with `--help/-h` option.

## Test run
To ensure that TOGA2 was installed properly, run the following command:

```bash
toga2.py test
```
Provided the successful installation, you will see the TOGA2 execution log printed to stdout, with the results stored at TOGA2/sample_output.
>[!NOTE]
> If you are running TOGA2 from Apptainer, you might have to provide a custom output directory with `--output/-o` option to bypass the read-only container configuration.

