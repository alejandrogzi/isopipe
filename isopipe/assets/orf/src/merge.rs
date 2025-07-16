use crate::{cli::MergeArgs, utils::*};
use config::{OverlapType, bed_to_struct_collection};
use packbed::{reader, record::GenePred};

use std::sync::Arc;

pub fn merge(args: MergeArgs) {
    // 1. need to read original bed file [DONE]
    // 2. read the orf and tai files -> use coordinates here to get CDS coordinates using bed file

    let records = bed_to_struct_collection::<GenePred>(
        reader(args.alignments)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));
}
