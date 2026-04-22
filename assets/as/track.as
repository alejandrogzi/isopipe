##########################################################
# Isopipe track description
##########################################################
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

##########################################################
# Isopipe bigWig track description [ spliceAi + aparent ]
##########################################################
track {BIGWIG_TRACK}
compositeTrack on
shortLabel HL ISOPIPE bigWig track
longLabel HL Isopipe bigWig track
group genes
priority 2
visibility full
itemRgb on
type bigWig
searchPriority 2.07207

        track {PREFIX}.aparent.forward
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL aparent forward
        longLabel HL_aparent_forward
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.aparent.forward.bw
        color 58,130,27

        track {PREFIX}.aparent.reverse
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL aparent reverse
        longLabel HL_aparent_reverse
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.aparent.reverse.bw
        color 212,154,25

        track {PREFIX}.spliceai.donor.forward
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL spliceai donor forward
        longLabel HL_spliceai_donor_forward
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.spliceai.donor.forward.bw
        color 38,53,201

        track {PREFIX}.spliceai.donor.reverse
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL spliceai donor reverse
        longLabel HL_spliceai_donor_reverse
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.spliceai.donor.reverse.bw
        color 114,33,148

        track {PREFIX}.spliceai.acceptor.forward 
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL spliceai acceptor forward
        longLabel HL_spliceai_acceptor_forward
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.spliceai.acceptor.forward.bw
        color 38,53,201

        track {PREFIX}.spliceai.acceptor.reverse
        parent {BIGWIG_TRACK}
        subtrack {BIGWIG_TRACK}
        shortLabel HL aparent forward
        longLabel HL_aparent_forward
        group genes
        priority 1
        visibility full
        itemRgb On
        type bigWig
        bigDataUrl {BROWSER}/{SPECIES}/{PREFIX}.aparent.forward.bw
        color 114,33,148
