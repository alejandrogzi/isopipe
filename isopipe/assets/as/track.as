# Isopipe track description
track HLIsoClassAnnot
compositeTrack on
shortLabel HL ISOPIPE annotations
longLabel Isopipe annotation track
group genes
priority 2
visibility pack
itemRgb on
type bigBed 12 +
searchPriority 2.07207
searchType bigBed

        track HLIsoClassAnnotPass
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL pass
        longLabel HL_Pass
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.pass.bb
        searchType bigBed

        track HLIsoClassAnnotFusions
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL fusions
        longLabel HL_Fusions
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.fusions.bb
        searchType bigBed

        track HLIsoClassAnnotRetention
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL retention
        longLabel HL_Retention
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.retention.bb
        searchType bigBed

        track HLIsoClassAnnotIntrprimming
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL intraprimming
        longLabel HL_Intraprimming
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.intraprimming.bb
        searchType bigBed

        track HLIsoClassAnnotNmd
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL nmd
        longLabel HL_Nmd
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.nmd.bb
        searchType bigBed

        track HLIsoClassAnnotTruncation
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL truncation
        longLabel HL_Truncation
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.pass.bb
        searchType bigBed

        track HLIsoClassAnnotRT
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL rt
        longLabel HL_rt
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.rt_reads.bb
        searchType bigBed

        track HLIsoClassAnnotTrash
        parent HLIsoClassAnnot
        subtrack HLIsoClassAnnot
        shortLabel HL trash
        longLabel HL_trash
        group genes
        priority 1
        visibility pack
        itemRgb On
        type bed 12 +
        bigDataUrl {BROWSER}/{SPECIES}/HLIsoClassAnnot.trash.bb
        searchType bigBed
