# Isopipe track description
track HLIsoClassAnnot
compositeTrack on
shortLabel HL ISOPIPE annotations
longLabel Isopipe annotation track
group genes
priority 2
visibility pack
itemRgb on
type bigBed 37
searchPriority 2.07207

        track HLIsoClassAnnotPass
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL pass
        longLabel HL_Pass
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.pass.bb

        track HLIsoClassAnnotFusions
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL fusions
        longLabel HL_Fusions
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.fusions.bb

        track HLIsoClassAnnotRetention
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL retention
        longLabel HL_Retention
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.retention.bb

        track HLIsoClassAnnotIntrprimming
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL intraprimming
        longLabel HL_Intraprimming
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.intraprimming.bb

        track HLIsoClassAnnotNmd
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL nmd
        longLabel HL_Nmd
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.nmd.bb

        track HLIsoClassAnnotTruncation
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL truncation
        longLabel HL_Truncation
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.truncation.bb

        track HLIsoClassAnnotRT
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL rt
        longLabel HL_rt
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.rt_reads.bb

        track HLIsoClassAnnotTrash
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL trash
        longLabel HL_trash
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bigBed 37
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.trash.bb
