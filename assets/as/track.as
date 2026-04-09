# Isopipe track description
track {TRACK}
compositeTrack on
shortLabel HL ISOPIPE annotations
longLabel HL Isopipe annotation track
group genes
priority 2
visibility pack
itemRgb on
type bigBed 37
searchPriority 2.07207

        track {PREFIX}.pass
        parent {TRACK}
        subtrack {TRACK}        
        shortLabel HL pass
        longLabel HL_Pass
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.pass.bb

        track {PREFIX}.retention
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL retention
        longLabel HL_Retention
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.retention.bb

        track {PREFIX}.intraprimming
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL intraprimming
        longLabel HL_Intraprimming
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.intraprimming.bb

        track {PREFIX}.truncation
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL truncation
        longLabel HL_Truncation
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.truncation.bb

        track {PREFIX}.rt
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL rt
        longLabel HL_rt
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.rt.bb

        track {PREFIX}.trash
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL trash
        longLabel HL_trash
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.trash.bb

        track {PREFIX}.orphans
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL orphans
        longLabel HL_orphans
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.orphans.bb

        track {PREFIX}.duplicates
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL duplicates
        longLabel HL_duplicates
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed {ADDITIONAL_COLUMNS}
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.duplicates.bb

        track {PREFIX}.fusions
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL fusions
        longLabel HL_fusions
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 12
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.fusions.bb

        track {PREFIX}.nmd
        parent {TRACK}
        subtrack {TRACK}
        shortLabel HL nmd
        longLabel HL_nmd
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 12
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.nmd.bb

