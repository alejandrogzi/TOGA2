use bigtools::BigBedRead;
use clap::Parser;
use fxhash::{FxHashMap, FxHashSet};
use std::cmp::{max, min};
use std::io::{BufRead, BufReader, Write, stdout};
use std::fs::File;
use std::sync::Arc;

use svgdom::{
    AttributeId, AttributeValue, Document, 
    ElementId, Length, LengthList, LengthUnit, 
    Node, NodeType, ViewBox
};
use svg::write;

const STDOUT: &str = "stdout";
const FREE_SPACE: f64 = 1.0;

struct Coords{
    start: u32,
    end: u32
}

struct Dimensions{
    height: f64,
    width: f64
}

fn recursive_addition<'a>(
    // from: &mut Document, 
    to: &mut Document,
    donor_node: &'a mut Node,
    recipient_node: &'a mut Node,
    dims: &mut Dimensions,
    mut height_mod: f64,
    plot_counter: u8
) {
    let mut is_svg_node: bool = false;
    for mut child in donor_node.children() {
        if child.node_type() == NodeType::Element {
            is_svg_node |= child.tag_id().unwrap() == ElementId::Svg;
            let mut new_node = to.create_element(
                child.tag_id().unwrap().clone()
            );
            // iterate over the elements
            for attr in  child.attributes().iter() {
                // ViewBox attribute; in TOGA2, projection-wise plots are guaranteed 
                // to have only one ViewBox 
                if let AttributeValue::ViewBox(view) = attr.value {
                    let box_height = view.h;
                    if plot_counter > 0 {dims.height += FREE_SPACE}
                    height_mod = dims.height;
                    dims.height += box_height;
                    let box_width = view.w;
                    dims.width = dims.width.max(box_width);
                    continue
                }
                if child.is_root() {
                    continue
                }
                if plot_counter > 1 {
                    let attr_id = attr.id().unwrap();
                    if attr_id == AttributeId::Y || attr_id == AttributeId::Y1 || attr_id == AttributeId::Y2 {
                        let adj_len = match &attr.value {
                            AttributeValue::Length(x) => {
                                AttributeValue::Length(Length::new(x.num + height_mod, LengthUnit::None))
                            },
                            AttributeValue::LengthList(x) => {
                                AttributeValue::LengthList(
                                    LengthList(
                                        x
                                        .iter()
                                        .map(|f| Length::new(f.num + height_mod, LengthUnit::None))
                                        .collect::<Vec::<Length>>()
                                    )
                                )
                            }
                            _ => {
                                continue
                            }
                        };
                        new_node.set_attribute((attr_id, adj_len));
                        continue
                    }
                }
                new_node.set_attribute(attr.clone());
            }
            if child.has_children() {
                recursive_addition(
                    to, 
                    &mut child,
                    if is_svg_node {recipient_node}  else {&mut new_node}, 
                    dims, 
                    height_mod, 
                    plot_counter
                );
            }
            if !is_svg_node {recipient_node.append(new_node.clone())};
        } else if child.node_type() == NodeType::Text {
            let new_node = to.create_node(NodeType::Text,child.text().clone());
            recipient_node.append(new_node.clone());
        }
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Extracts SVG plots for given projections from a TOGA2 BigBed file
struct Args{
    /// Input BigBed file
    #[arg()]
    bigbed_file: String,

    /// Bed file corresponding to the input BigBed
    #[arg()]
    bed_file: String,

    /// A comma-separated list of projections to fetch the plots for
    #[arg()]
    projections: String,

    /// If set, treats input names as transcript names
    #[arg(long, short = 't', default_value_t = false)]
    transcripts: bool,

    /// Path to output file; if no argument provided or if set to "stdout", the results 
    /// are printed to standard output
    #[arg(long, short = 'o', default_value_t = String::from(STDOUT))]
    output: String
}

fn main() {
    let args = Args::parse();
    // parse the input projections
    let projections: Vec<&str> = args.projections.split(',').collect();
    if projections.len() == 0 {panic!("No valid projection names provided")}
    let mut projections_found: FxHashSet<String> = FxHashSet::default();

    // open the guiding Bed file and extract the coordinates for each projection]
    let mut chrom2choords: FxHashMap<String, Coords> = FxHashMap::default();
    let bed_handle = {
        let path = File::open(args.bed_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for (i, line_) in bed_handle.lines().enumerate() {
        if let Ok(line) = line_ {
            // let bed_record = parse_bed(line, 12, true);
            let proj_comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if proj_comps.len() < 4 {
                panic!(
                    "Improper Bed file formatting at line {}; expected at least 4 columns, got {}",
                    i+1, proj_comps.len()
                )
            }
            let proj_name: &str = proj_comps[3];
            if !projections.contains(&proj_name) {continue}
            let chrom: String = String::from(proj_comps[0]);
            let start: u32 = proj_comps[1]
                .parse::<u32>()
                .expect(
                    &format!(
                        "Improper Bed file formatting at line {}: {} is not a valid 64-bit integer", 
                        i+1, proj_comps[1]
                    )
                );
            let end: u32 = proj_comps[2]
                .parse::<u32>()
                .expect(
                    &format!(
                        "Improper Bed file formatting at line {}: {} is not a valid 64-bit integer", 
                        i+1, proj_comps[2]
                    )
                );
            chrom2choords
                .entry(chrom)
                .and_modify(|x| {
                    x.start = min(x.start, start);
                    x.end = max(x.end, end)
                })
                .or_insert(Coords{start, end});
        }
    }

    // create a combined plot stub
    let mut combined_plot = Document::new();
    let mut svg = combined_plot.create_element(ElementId::Svg);
    combined_plot.root().append(svg.clone());
    let mut dims = Dimensions { height: 0.0, width: 0.0 };
    let mut plot_counter: u8 = 1;
    let mut height_mod: f64 = 0.0;
    // for each requested projection, fetch its entry from the BigBed file
    if let Ok(mut bb) = BigBedRead::open_file(args.bigbed_file) {
        for (chrom, coords) in chrom2choords.iter() {
            if let Ok(entry_inter) = bb.get_interval(chrom, coords.start, coords.end) {
                for entry_ in entry_inter {
                    if let Ok(entry) = entry_ {
                        let rest_comps = entry.rest.split('\t').collect::<Vec<&str>>();
                        let name: &str = rest_comps[0];
                        if !projections.contains(&name) {continue}
                        // save the plot string
                        let plot = rest_comps[27];
                        // parse the SVG string as a document
                        let tree = Document::from_str(plot).unwrap();
                        recursive_addition(
                            &mut combined_plot, 
                            &mut tree.root(),//svg_element().unwrap(),//tree.root(), 
                            &mut svg, 
                            &mut dims, 
                            height_mod, 
                            plot_counter
                        );
                        projections_found.insert(String::from(name));
                        plot_counter += 1;
                    }
                }
            }
        }
    }

    // println!("Time to plot!");

    svg.set_attribute(
        (
            AttributeId::ViewBox, 
            AttributeValue::ViewBox(ViewBox::new(0.0, 0.0, dims.width, dims.height))
        )
    );
    svg.set_attribute(
        (
            AttributeId::Height,
            AttributeValue::Length(Length::new(dims.height, LengthUnit::None)),
        )
    );
    svg.set_attribute(
        (
            AttributeId::Width,
            AttributeValue::Length(Length::new(dims.width, LengthUnit::None)),
        )
    );
    for child in combined_plot.root().descendants() {
        // println!("root child: {:#?}", child);
        for attr in child.attributes().iter() {
            // println!("attr={:#?}", attr);
        }
    }
    // println!("The first SVG-tagged element: {:#?}", combined_plot.svg_element());
    println!("Combined plot: \n{}\n", combined_plot.to_string());
}