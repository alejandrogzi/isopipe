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
    ret->R_read_intron_status = cloneString(row[2]);
    ret->R_metadata_html = cloneString(row[3]);
    ret->R_retentions_html = cloneString(row[4]);
    ret->R_retains_rt_intron_html = cloneString(row[5]);
    ret->R_rt_html = cloneString(row[6]);
    ret->P_read_status = cloneString(row[7]);
    ret->P_metadata_html = cloneString(row[8]);
    ret->T_read_status = cloneString(row[9]);
    ret->T_metadata_html = cloneString(row[10]);
    ret->O_read_orf_score = cloneString(row[11]);
    ret->O_metadata_html = cloneString(row[12]);
    ret->TAG_five_clip_len = cloneString(row[13]);
    ret->TAG_three_clip_len = cloneString(row[14]);
    ret->TAG_polya_len = cloneString(row[15]);
    ret->TAG_polya_read_len = cloneString(row[16]);
    ret->TAG_mapping_identity = cloneString(row[17]);
    ret->TAG_singleton = cloneString(row[18]);
    ret->TAG_orf_number = cloneString(row[19]);
    ret->TAG_nested_orf_number = cloneString(row[20]);
    ret->TAG_fake_fusion = cloneString(row[21]);
    ret->TAG_review_fusion = cloneString(row[22]);
    ret->TAG_strong_nmd = cloneString(row[23]);
    ret->TAG_weak_nmd = cloneString(row[24]);
    ret->TAG_unique_tai = cloneString(row[25]);
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
       "  <li>Has RT intron: %s</li>\n"
       "  <li>Intrapriming status: %s</li>\n"
       "  <li>Truncation status: %s</li>\n"
       "  <li>ORF prediction score: %s</li>\n"
       "</ul><br>\n",
       info->R_read_status,
       info->R_read_intron_status,
       info->P_read_status,
       info->T_read_status,
       info->O_read_orf_score);

htmlHorizontalLine();

//printf("<BR><a data-toggle=\"collapse\" href=\"#collapseIR\">Show details of intron retention detection</a>\n");
//printf("<div id=\"collapseIR\" class=\"panel-collapse collapse\">\n");

// Tag details

printf("<details open>\n"
       "  <summary>Show read tag details</summary>\n"
       "  <ul>\n"
       "    <li>TAG_five_clip_len: %s</li>\n"
       "    <li>TAG_three_clip_len: %s</li>\n"
       "    <li>TAG_polya_len: %s</li>\n"
       "    <li>TAG_polya_read_len: %s</li>\n"
       "    <li>TAG_mapping_identity: %s</li>\n"
       "    <li>TAG_singleton: %s</li>\n"
       "    <li>TAG_orf_number: %s</li>\n"
       "    <li>TAG_nested_orf_number: %s</li>\n"
       "    <li>TAG_fake_fusion: %s</li>\n"
       "    <li>TAG_review_fusion: %s</li>\n"
       "    <li>TAG_unique_tai: %s</li>\n"
       "  </ul>\n"
       "</details>\n",
       info->TAG_five_clip_len,
       info->TAG_three_clip_len,
       info->TAG_polya_len,
       info->TAG_polya_read_len,
       info->TAG_mapping_identity,
       info->TAG_singleton,
       info->TAG_orf_number,
       info->TAG_nested_orf_number,
       info->TAG_fake_fusion,
       info->TAG_review_fusion,
       info->TAG_unique_tai);

htmlHorizontalLine();

// Intron retention details

printf("<details open>\n"
       "  <summary>Show details of read intron retention</summary>\n"
       "  <h4>Read status details</h4>\n"
       "  %s\n"
       "\n"
       "  <h4>Retention details</h4>\n"
       "  %s\n"
       "\n"
       "  <h4>Retained RT intron details</h4>\n"
       "  %s\n"
       "\n"
       "  <h4>Spliced RT intron details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->R_metadata_html,
       info->R_retentions_html,
       info->R_retains_rt_intron_html,
       info->R_has_rt_html);

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
       info->P_metadata_html);

htmlHorizontalLine();

// ORF prediction details

printf("<details open>\n"
       "  <summary>Show details of ORf prediction</summary>\n"
       "  <h4>ORF details</h4>\n"
       "  %s\n"
       "</details>\n",
       info->O_metadata_html);


hPrintf("<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css\">");
hPrintf("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js\"></script>");
hPrintf("<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js\"></script>");

printTrackHtml(tdb);
}
