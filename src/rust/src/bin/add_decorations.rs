use clap::Parser;
use fxhash::FxHashMap;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Given a Bed file and a TOGA2 mutation report file, prepare 
/// 
struct Args {

    /// path to the Bed12 file to decorate
    #[arg(long, short = 'b')]
    bed_file: String,

    ///  path to mutation file
    #[arg(long, short = 'm')]
    mutation_file: String,

    /// path to the output file
    #[arg(long, short = 'o')]
    output: String
}

// static MUT2SHAPE: FxHashMap<&str, String> = FxHashMap::from_iter([
//     ("", String::from("")),
// ]);

fn main() {
    let args = Args::parse();
}