use bed2gtf;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(version = "1.0", about, long_about = None)]
/// Converts file
struct Args {
    /// Input Bed file
    #[arg(long, short = 'i', default_value_t = String::from("stdout"))]
    input: String
}

fn main() {}