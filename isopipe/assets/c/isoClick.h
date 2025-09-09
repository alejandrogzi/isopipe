/* isoClassClick - click handling for IsoSeq classification tracks */

/*
 * author: Alejandro Gonzales-Irribarren, 2025
 * email = "alejandrxgzi@gmail.com"
 * github = "https://github.com/alejandrogzi"
 * version = "0.0.1"
 */

#ifndef ISOCLASSCLICK_H
#define ISOCLASSCLICK_H

#define YES_ "Yes"
#define NO_ "No"
#define ONE_ "1"

#define HLISO_BED_PREFIX "HLIsoClass"
#define HLISO_BED_PREFIX_LEN 11

#define HLISO_MAXCHAR 255

struct isoClassDataBB
{
    char* R_read_status;                // R: Read status
    char* R_read_intron_status;         // R: Intron status
    char* R_metadata_html;              // R: Metadata in HTML
    char* R_retentions_html;            // R: Retention introns in HTML
    char* R_retains_rt_intron_html;     // R: RT intron retention in HTML
    char* R_rt_html;                    // R: Reverse transcription events in HTML
    char* P_read_status;                // P: Read status
    char* P_metadata_html;              // P: Metadata in HTML
    char* T_read_status;                // T: Read status
    char* T_metadata_html;              // T: Metadata in HTML
    char* O_read_orf_score;             // O: ORF score
    char* O_metadata_html;              // O: Metadata in HTML
    char* TAG_five_clip_len;            // Length of 5' clipping
    char* TAG_three_clip_len;           // Length of 3' clipping
    char* TAG_polya_len;                // PolyA tail length
    char* TAG_polya_read_len;           // Length of read contributing to PolyA
    char* TAG_mapping_identity;         // Mapping identity (%)
    char* TAG_singleton;                // Singleton read flag
    char* TAG_orf_number;               // Number of ORFs
    char* TAG_nested_orf_number;        // Number of nested ORFs
    char* TAG_fake_fusion;              // Fake fusion flag
    char* TAG_review_fusion;            // Fusion under review flag
    char* TAG_strong_nmd;               // Strong NMD prediction
    char* TAG_weak_nmd;                 // Weak NMD prediction
    char* TAG_unique_tai;               // Unique translationAi prediction
};


struct isoClassDataBB *isoClassDataBBLoad(char **row, bits16 fieldCount);
/* Load a isoClassData from bigBed. Dispose of this with isoClassDataBBFree(). */

void isoClassDataBBFree(struct isoClassDataBB **pEl);
/* Free a single dynamically allocated isoClassData such as created
 * with isoClassDataBBLoad. */

void doHillerLabIsoClass(char *database, struct trackDb *tdb, char *item, char *table_name);
/* Put up IsoSeq classification track info. */

#endif  // ISOCLASSCLICK_H
