table isopipe
"Schema for base isopipe-produced long-read transcriptomic reads [fusions + nmd]"
(
    string  chrom;              "Reference sequence chromosome or scaffold"
    uint    chromStart;         "Start position of feature on chromosome"
    uint    chromEnd;           "End position of feature on chromosome"
    string  name;               "Feature ID"
    uint    score;              "Score (0-1000)"
    char[1] strand;             "+ or - for strand"
    uint    thickStart;         "Coding region start"
    uint    thickEnd;           "Coding region end"
    uint    itemRgb;           "RGB color value"
    int     blockCount;         "Number of exons"
    int[blockCount] blockSizes; "Exon sizes"
    int[blockCount] chromStarts;"Exon starts relative to chromStart"
)
