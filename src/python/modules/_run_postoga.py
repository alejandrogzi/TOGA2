#!/beegfs/projects/hillerlab/genome/bin/scripts/postoga/.env/bin/python3

import argparse
import os
import shutil
import subprocess
from typing import Union

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

REFERENCES = ["hg38", "mm10", "galGal6", "HLbosTau10", "HLeleMaxInd3A"]
FORMATS = ["gtf", "gff", "bed"]
POSTOGA = "/beegfs/projects/hillerlab/genome/bin/scripts/postoga/postoga.py"
TOML = f"{POSTOGA.replace('postoga.py', 'pyproject.toml')}"
POSTOGA = "/beegfs/projects/hillerlab/genome/bin/scripts/postoga/postoga.py"
VENV = "/beegfs/projects/hillerlab/genome/bin/scripts/postoga/.venv/bin/activate"


def main() -> None:
    __hi()
    __check()

    args = parse_args()
    run(args)


def __get_query_directory(
    ref: Union[str, os.PathLike[str]], query: Union[str, os.PathLike[str]]
) -> Union[str, os.PathLike[str]]:
    """
    Get TOGA directory for given assemblies.

    Parameters
    ----------
    ref : str
        Reference assembly ID.
    query : str
        Query assembly ID.

    Returns
    -------
    str
        Path to TOGA directory

    Examples
    --------
    >>> _get_toga_dir_for_asms("hg38", "dasNov3")
    >>> "/projects/hillerlab/genome/gbdb-HL/hg38/TOGA2/vs_dasNov3"
    """
    return f"/projects/hillerlab/genome/gbdb-HL/{ref}/TOGA2/vs_{query}"


def run(args: argparse.Namespace) -> None:
    """
    Run postoga on TOGA2 results.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments

    Returns
    -------
    None

    Examples
    --------
    >>> run(args)
    """
    if args.toga_dir:
        query_dir = args.toga_dir
    else:
        query_dir = __get_query_directory(args.reference, args.query)

    postoga_dir = os.path.join(query_dir, args.dir)

    # INFO: if overwrite, delete postoga subdir + .gtf
    if args.force:
        try:
            shutil.rmtree(postoga_dir)
            os.remove(os.path.join(query_dir, "query_annotation.gtf"))
        except FileNotFoundError:
            print(
                "WARN: no postoga subdir nor query_annotation.gtf found, continuing..."
            )

    if args.env:
        cmd = f"source {args.env} && {POSTOGA} base --togadir {query_dir} --outdir {postoga_dir} --to {args.to} --target {args.target}"
    else:
        cmd = f"source {VENV} && {POSTOGA} base --togadir {query_dir} --outdir {postoga_dir} --to {args.to} --target {args.target}"

    print(f"INFO: Running {cmd}")

    sh = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        check=True,
        cwd="/beegfs/projects/hillerlab/genome/bin/scripts/postoga",
    )
    print(sh.stdout)

    # INFO: move postoga results to the TOGA directory
    postoga_files = os.listdir(postoga_dir)
    for file in postoga_files:
        if file.endswith("gtf.gz") or file.endswith("gff.gz") or file.endswith("table"):
            print(f"INFO: moving {file}")

            if file.endswith("gtf.gz"):
                shutil.move(os.path.join(postoga_dir, file), os.path.join(query_dir, "query_annotation.gtf.gz"))
            else:
                shutil.move(os.path.join(postoga_dir, file), query_dir)

    print(f"INFO: removing {postoga_dir}")
    shutil.rmtree(postoga_dir, ignore_errors=True)

    return


def __hi() -> None:
    """
    Print a welcome message.

    Returns
    -------
    None

    Examples
    --------
    >>> __hi()
    """
    start_message = "\n> toga2_run_postoga.py"
    version = f"> version: {__version__}\n\n"

    print(start_message + "\n" + version + "\n\n")

    return


def __check() -> None:
    """
    Check if all required programs are installed.

    Returns
    -------
    None

    Examples
    --------
    >>> check()
    """
    if not os.path.exists(POSTOGA):
        raise FileNotFoundError(f"ERROR: {POSTOGA} not found!")

    if not os.path.exists(TOML):
        raise FileNotFoundError("ERROR: pyproject.toml not found!")

    print("INFO: All required programs are installed!")

    return


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments

    Examples
    --------
    >>> parse_args()
    Namespace(query=None, ref=None, version='0.0.1')
    """

    app = argparse.ArgumentParser("Run postoga on TOGA2 results")

    app.add_argument(
        "-r",
        "--reference",
        help=f"Reference assembly ID; options: {REFERENCES}",
        required=True,
        type=str,
        choices=REFERENCES,
    )
    app.add_argument(
        "-q",
        "--query",
        help="Query assembly ID, for example: dasNov3",
        required=True,
        type=str,
    )
    app.add_argument(
        "-td",
        "--toga_dir",
        help="Overrides Hillerlab common path to TOGA2 dirs and specify a custom dir",
        required=False,
        type=str,
    )
    app.add_argument(
        "-v", "--version", required=False, default=f"{__version__}", help="Version"
    )
    app.add_argument(
        "-t",
        "--to",
        help="Format to convert TOGA .bed files to; options: gtf, gff, bed",
        choices=FORMATS,
        type=str,
        default="gtf",
    )
    app.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force overwriting",
    )
    app.add_argument(
        "-e",
        "--env",
        help="Virtual environment with packages [activate file!]",
        type=str,
    )
    app.add_argument(
        "-d",
        "--dir",
        help="Overrides TARGET/postoga as the new output postoga dir",
        type=str,
        default="postoga",
    )
    app.add_argument(
        "-tg",
        "--target",
        help="Specify the .bed input file to used by the program",
        type=str,
        choices=["bed", "utr", "both"],
        default="utr",
    )

    args = app.parse_args()

    return args


if __name__ == "__main__":
    main()
