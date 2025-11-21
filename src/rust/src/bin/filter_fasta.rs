use clap::Parser;
use fxhash::{FxHashSet};
use rust::read;
use std::io::{BufRead, stdout, Write};
use std::fs::File;
use std::path::Path;


#[derive(Parser, Debug)]
#[command(version = "1.0", about, long_about = None)]
/// Fasta filter for TOGA2; Given the pairwise alignment file
/// and the list of deprecated projections as a single-column list file, 
/// produces a file of trusted query sequences
struct Args {
    /// Input Fasta file; if no value is provided or if value set to "stdin", 
    /// reads the data from the standard input stream
    #[arg(long, short = 'i', default_value_t = String::from("stdin"))]
    input: String,

    /// A list of deprecated projections; 
    /// those will be filtered out from the final output
    #[arg(long, short = 'd')]
    deprecated: Option<String>,

    /// A final Bed file to filter the Fasta file by
    #[arg(long, short = 'b')]
    bed_file: Option<String>,

    /// A path to save the results to; if no value provided or if set to "stdout",
    /// the results are printed to standard output stream
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String
}

fn main() {
    let args = Args::parse();

    // if a list of deprecated projections was provided, copy its contents to a set
    // let mut excluded_projs: FxHashSet<String> = FxHashSet::default();
    // if let Some(deprecated) = args.deprecated {
    //     for line_ in read(deprecated).lines() {
    //         if let Ok(line) = line_ {
    //             excluded_projs.insert(line);
    //         }
    //     }
    // }
    let mut included_projs: FxHashSet<String> = FxHashSet::default();
    let mut paralogs: FxHashSet<String> = FxHashSet::default();
    let mut ppgenes: FxHashSet<String> = FxHashSet::default();
    if let Some(bed_file) = args.bed_file {
        for line_ in read(bed_file).lines() {
            if let Ok(line) = line_ {
                let comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
                let name = comps[3]
                    .split('$') // for fragmented projections, ignore the the fragment number
                    .collect::<Vec<&str>>()
                    .first()
                    .unwrap()
                    .replace("#retro", "")// remove the retrogene postfix
                    .replace("#paralog", ""); // remove the paralog postfix
                included_projs.insert(name.clone());
                if comps[3].contains("#paralog") {paralogs.insert(name.clone());}
                if comps[3].contains("#retro") {ppgenes.insert(name.clone());}
            }
        }
    }

    // create an output handle
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    // read and filter the input Fasta
    let mut header: String = String::new();
    let mut seq: String = String::new();
    let mut record: bool = false;
    for line_ in read(args.input).lines() {
        if let Ok(line) = line_ {
            if line.starts_with('>') {
                // a full Fasta entry has been recorded; write it to the output file
                if record {
                    if let Err(e) = output_file.write(format!("{}\n{}\n", header, seq).as_bytes()) {
                        eprintln!("Failed to write the line: {}", e);
                    }
                    header = String::new();
                    seq = String::new();
                }
                let comps: Vec<String> = line
                    .split('|')
                    .map(|x| x.replace(' ', ""))
                    .collect::<Vec<String>>();
                // consistent with the TOGA2 Fasta header notation, we assume that all the sequences 
                // either represent query by default or have their source specified in the last header field
                if comps.len() > 1 {
                    let is_query: bool = comps[comps.len() - 1] == "QUERY";
                    // do not record the sequence line if it does not belong to query
                    if !is_query {record = false; continue}
                }
                let proj = comps[0]
                    .replace('>', "")
                    .replace("#retro", "")
                    .replace("#paralog", "");
                // do not record the deprecated projection's sequence
                // if excluded_projs.contains(&proj) {record = false; continue}
                if included_projs.len() > 0 && !included_projs.contains(&proj) {continue}
                header = comps[0].clone();
                if paralogs.contains(&proj) {header += "#paralog"}
                if ppgenes.contains(&proj) {header += "#retro"}
                record = true;
                // proceed to the next line
                continue
            }
            if record {
                seq += &line.replace('-', "");
            }
        }
    }
    // write the final entry to the file
    if record {
        if let Err(e) = output_file.write(format!("{}\n{}", header, seq).as_bytes()) {
            eprintln!("Failed to write the line: {}", e);
        }
    }
}