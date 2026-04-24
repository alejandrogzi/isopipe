/* isoClick - click handling for IsoSeq classification tracks */

/*
 * author: Alejandro Gonzales-Irribarren, 2025
 * credits: Michael Hiller, Bogdan Kirilenko
 * email = "alejandrxgzi@gmail.com"
 * github = "https://github.com/alejandrogzi"
 * version = "0.0.1"
*/

#include "common.h"
#include "hgc.h"
#include "isoClick.h"
#include "string.h"
#include "htmshell.h"
#include "chromAlias.h"


struct isoClassDataBB *isoClassDataBBLoad(char **row, bits16 fieldCount)
/* Load a isoClassData from bigBed.  Dispose of this with isoClassDataFree(). */
{
    struct isoClassDataBB *ret;
    AllocVar(ret);
    ret->R_read_status = cloneString(row[1]);
    ret->R_read_code = cloneString(row[2]);
    ret->R_metadata_html = cloneString(row[3]);
    ret->T_read_status = cloneString(row[4]);
    ret->T_read_code = cloneString(row[5]);
    ret->T_metadata_html = cloneString(row[6]);
    ret->I_read_status = cloneString(row[7]);
    ret->I_read_code = cloneString(row[8]);
    ret->I_metadata_html = cloneString(row[9]);
    ret->O_metadata_html = cloneString(row[10]);
    ret->C_collapsed = cloneString(row[11]);

    return ret;
}


void isoClassDataBBFree(struct isoClassDataBB **pEl)
/* Free a single dynamically allocated isoClassDatasuch as created
 * with isoClassDataLoad(). */
{}

void doHillerLabIsoClass(char *database, struct trackDb *tdb, char *item, char *table_name)
/* Put up isoClass Gene track info. */
{
int start = cartInt(cart, "o");
int end = cartInt(cart, "t");
char *chrom = cartString(cart, "c");

/* skip if not a bigBed */
if (! (startsWith("bigBed", tdb->type)))
        return;

char *fileName = bbiNameFromSettingOrTable(tdb, NULL, tdb->table);
struct bbiFile *bbi =  bigBedFileOpenAlias(hReplaceGbdb(fileName), chromAliasFindAliases);
struct lm *lm = lmInit(0);
struct bigBedInterval *bbList = bigBedIntervalQuery(bbi, chrom, start, end, 0, lm);
struct bigBedInterval *bb;
char *fields[bbi->fieldCount];
for (bb = bbList; bb != NULL; bb = bb->next)
    {
    if (!(bb->start == start && bb->end == end))
        continue;

    char *name = cloneFirstWordByDelimiterNoSkip(bb->rest, '\t');
    boolean match = (isEmpty(name) && isEmpty(item)) || sameOk(name, item);
    if (!match)
        continue;

    char startBuf[16], endBuf[16];
    bigBedIntervalToRow(bb, chrom, startBuf, endBuf, fields, bbi->fieldCount);
    break;
    }

printf("<h3>Read: %s</h3>\n", item);
struct isoClassDataBB *info = isoClassDataBBLoad(&fields[11], bbi->fieldCount);  // Bogdan: why 11? 0-11 are bed-like fields likely

printf("<h4>Overall classification</h4>\n"
       "<ul>\n"
       "  <li>Intron retention status: %s</li>\n"
       "  <li>Intrapriming status: %s</li>\n"
       "  <li>Truncation status: %s</li>\n"
       "</ul><br>\n",
       info->R_read_status,
       info->I_read_status,
       info->T_read_status
);

htmlHorizontalLine();

//printf("<BR><a data-toggle=\"collapse\" href=\"#collapseIR\">Show details of intron retention detection</a>\n");
//printf("<div id=\"collapseIR\" class=\"panel-collapse collapse\">\n");

// Intron retention details

printf("<details open>\n"
       "  <summary>Show details of read intron retention</summary>\n"
       "  <h4>Read status details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->R_metadata_html);

htmlHorizontalLine();

// Truncation details

printf("<details open>\n"
       "  <summary>Show details of read start truncation</summary>\n"
       "  <h4>Truncation details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->T_metadata_html);

htmlHorizontalLine();

// Intraprimming details

printf("<details open>\n"
       "  <summary>Show details of read end intraprimming</summary>\n"
       "  <h4>Intraprimming details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->I_metadata_html);

htmlHorizontalLine();

// ORF prediction details

printf("<details open>\n"
       "  <summary>Show details of ORF prediction</summary>\n"
       "  <h4>ORF details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->O_metadata_html);

htmlHorizontalLine();

// Collapsed details

printf("<details open>\n"
       "  <summary>Show details of collapse queue</summary>\n"
       "  <h4>Collapse queue details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->C_collapsed);


hPrintf("<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css\">");
hPrintf("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js\"></script>");
hPrintf("<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js\"></script>");

printTrackHtml(tdb);
}
