/* isoClick - click handling for IsoSeq classification tracks */

/*
 * author: Alejandro Gonzales-Irribarren, 2025
 * credits: Michael Hiller, Bogdan Kirilenko
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
#define HLISO_BED_PREFIX_LEN 10

#define HLISO_MAXCHAR 255

struct isoClassDataBB
{
    char* R_read_status;                // R: Read status
    char* R_read_code;                  // R: Read code
    char* R_metadata_html;              // R: Metadata in HTML
    char* T_read_status;                // T: Read status
    char* T_read_code;                  // T: Read code
    char* T_metadata_html;              // T: Metadata in HTML
    char* I_read_status;                // P: Read status
    char* I_read_code;                  // P: Read code
    char* I_metadata_html;              // P: Metadata in HTML
    char* O_metadata_html;              // O: Metadata in HTML
    char* C_collapsed;                  // C: Collapsed list of reads
};


struct isoClassDataBB *isoClassDataBBLoad(char **row, bits16 fieldCount);
/* Load a isoClassData from bigBed. Dispose of this with isoClassDataBBFree(). */

void isoClassDataBBFree(struct isoClassDataBB **pEl);
/* Free a single dynamically allocated isoClassData such as created
 * with isoClassDataBBLoad. */

void doHillerLabIsoClass(char *database, struct trackDb *tdb, char *item, char *table_name);
/* Put up IsoSeq classification track info. */

#endif  // ISOCLASSCLICK_H
