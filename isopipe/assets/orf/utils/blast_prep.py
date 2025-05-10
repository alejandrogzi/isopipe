"""
Prepare data and run BLASTP & TransDecoder
"""

import logging
import shutil
import subprocess
import copy
import os

from .common import chunk_fasta, multiplex_dataset
from bedutils.bed import BedRow
from .merging import reformat_id_orfipy

logger = logging.getLogger(__name__)

NF_FILE_LINK = (
    "/beegfs/projects/hillerlab/genome/src/ORFTree/nextflow/chunked_diamond.nf"
)
NF_CONFIG = "/beegfs/projects/hillerlab/genome/src/ORFTree/nextflow/nextflow.config"


def feed_toga(toga_bed, assembly, tmp_dir):
    with open(toga_bed, "r") as f:
        # Convert to BED6
        with open(f"{tmp_dir}/ORF/TOGA_BED.bed", "w") as bed_output:
            for line in f:
                cols = line.split("\t")
                cols[3] = f"TOGA_cols{cols[3]}"
                bed_output.write("\t".join(cols[:5]) + "\n")

    # Extract raw sequences for TOGA
    subprocess.run(
        f"twoBitToFa {assembly} {tmp_dir}/ORF/TOGA_READS.fa -bed={tmp_dir}/ORF/TOGA_BED.bed"
    )


def de_multiplex_blast(index_path, blast_results, output_path):
    index_dict = dict()
    with open(index_path, "r", encoding="utf-8") as index_file:
        for line in index_file:
            cols = line.strip().split("\t")
            ids = cols[1].split("⍼")
            index_dict[cols[0]] = ids

    with open(blast_results, "r") as blast_file:
        with open(output_path, "w") as output_file:
            for line in blast_file:
                # Note: no strip() needed here
                cols = line.split("\t")
                multiplexed_id = cols[0]
                try:
                    ids = index_dict.get(multiplexed_id)

                    for my_id in ids:
                        cols[0] = my_id
                        output_file.write("\t".join(cols))

                except KeyError:
                    logger.warning(f"Key {multiplexed_id} not found in BLAST index")


def run_transdecoder_long_orfs(input_file, output_dir) -> bool:
    """
    Wrapper for TransDecoder Long ORF
    :param input_file: input file
    :param output_dir: output file
    :return: True if successful, False otherwise
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    transdecoder_args = [
        "TransDecoder.LongOrfs",
        "--output_dir",
        output_dir,
        "-m",
        "30",
        "-t",
        input_file,
    ]
    transdecoder_process = subprocess.run(transdecoder_args)

    if transdecoder_process.returncode != 0:
        raise RuntimeError(transdecoder_args)

    return True


def run_orfipy(input_file, output_dir) -> bool:
    """
    Wrapper for OrfiPy ORF
    :param input_file:
    :param output_dir:
    :return:
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    orf_args = [
        "orfipy",
        input_file,
        "--pep",
        "orfs_raw.pep",
        "--bed",
        "orfs.bed",
        "--partial-5",
        "--partial-3",
        "--include-stop",
        "--min",
        "100",
        "--ignore-case",
        "--outdir",
        output_dir,
    ]
    orf_process = subprocess.run(orf_args)

    if orf_process.returncode != 0:
        raise RuntimeError(orf_args)

    return True


def run_diamond(
    input_file,
    output_file,
    database,
    parent_tmp,
    num_threads=-1,
    chunk_size=-1,
    nextflow=False,
):
    """
    Run Diamond via para.
    :param parent_tmp: Parent temporary directory
    :param input_file:
    :param output_file:
    :param database:
    :param num_threads:
    :param chunk_size:
    :return:
    """

    if nextflow:
        command = f"nextflow {NF_FILE_LINK} -c {NF_CONFIG} --query {input_file} --db {database} --out {output_file} -work-dir {parent_tmp}/TEMP_DIAMOND"
        job_run = subprocess.run(command, shell=True)
        return job_run
    else:
        jobs = []
        chunked_tmp = f"{parent_tmp}/BLASTP"
        os.makedirs(f"{parent_tmp}/BLASTP")
        chunk_count = chunk_fasta(input_file, chunk_size, chunked_tmp, "BLASTP")

        logger.info(
            f"Diamond will run on {chunk_count} chunks with {chunk_size} sequences each."
        )

        # Create BLAST JOBS for chunked runs
        for i in range(chunk_count):
            jobs.append(
                f"diamond blastp -e 1e-10 --sensitive --query {chunked_tmp}/BLASTP_{i}.pep --db {database} -o {chunked_tmp}/BLASTP_{i}.txt\n"
            )

        with open(f"{chunked_tmp}/JOBFILE", "w") as jobfile:
            jobfile.writelines(jobs)

        # Run on cluster
        job_run = subprocess.run(
            f"para make BLASTP_{parent_tmp} {chunked_tmp}/JOBFILE -q shortmed -memoryMb 15000",
            shell=True,
            check=True,
        )

        # Collect results
        logger.info("Merging results from jobs into a single file...")
        subprocess.run(
            f"cat {chunked_tmp}/*.txt >> {output_file}", shell=True, check=True
        )

        return job_run


def run_blastp(
    input_file, output_file, database, parent_tmp, num_threads=16, chunk_size=-1
):
    """
    Run BLASTP via para.
    :param parent_tmp: Parent temporary directory
    :param input_file:
    :param output_file:
    :param database:
    :param num_threads:
    :param chunk_size:
    :return:
    """
    jobs = []
    chunked_tmp = f"{parent_tmp}/BLASTP"
    os.makedirs(f"{parent_tmp}/BLASTP")
    chunk_count = chunk_fasta(input_file, chunk_size, chunked_tmp, "BLASTP")

    logger.info(
        f"BLASTP will run on {chunk_count} chunks with {chunk_size} sequences each."
    )

    # Create BLAST JOBS for chunked runs
    for i in range(chunk_count):
        jobs.append(
            f"blastp -evalue 1e-10 -outfmt 6 -num_threads {num_threads} -db {database} -query {chunked_tmp}/BLASTP_{i}.pep -out {chunked_tmp}/BLASTP_{i}.txt\n"
        )

    with open(f"{chunked_tmp}/JOBFILE", "w") as jobfile:
        jobfile.writelines(jobs)

    # Run on cluster
    job_run = subprocess.run(
        f"para make BLASTP_{parent_tmp} {chunked_tmp}/JOBFILE -q shortmed",
        shell=True,
        check=True,
    )

    # Collect results
    logger.info("Merging results from jobs into a single file...")
    subprocess.run(f"cat {chunked_tmp}/*.txt >> {output_file}", shell=True, check=True)

    return job_run


def filter_inphase(in_file, out_file):
    with open(out_file, "w") as o:
        with open(in_file, "r") as i:
            for line in i:
                cols = line.strip().split("\t")
                if cols[-1] == "+":
                    o.write(line)


def export_cds_orfs_bed(
    bed_file, aligned_bed, nested_data, output_bed, tmp_dir, transdecoder=False
):
    """
    Export CDS and nested ORFs from a BED-based output of an ORF finder -> specifically orfipy
    :param bed_file: BED ORF file -> chroms are transcript IDs
    :param aligned_bed: The transcript alignmend BED file
    :param nested_data: .tsv nested ORF file
    :param output_bed: Output path
    :param tmp_dir: Unified directory for temporary files to use
    :param transdecoder: is the BED file transdecoder formatted
    :return:
    """

    cds_dict = dict()
    logger.info(f"Loading predicted CDS from {bed_file}")

    with open(bed_file, "r") as bed_file:
        for line in bed_file:
            cols = line.strip().split("\t")
            if transdecoder:
                if cols[7] != "CDS":
                    continue
                my_id = cols[3][4:]
                cds_dict[my_id] = {
                    "start": int(cols[1]),
                    "stop": int(cols[2]),
                    "strand": cols[6],
                }
            else:
                my_id = cols[3].split(";")[0][3:].replace("_ORF.", ".p")
                cds_dict[my_id] = {
                    "start": int(cols[1]),
                    "stop": int(cols[2]),
                    "strand": cols[5],
                }

    export_rows = []
    export_dict = dict()
    dropped_count = 0
    logger.info(f"Trimming {aligned_bed} to CDS coordinates")
    print(list(cds_dict.keys())[:10])
    missed = []
    with open(aligned_bed, "r") as ref_file:
        for line in ref_file:
            x = 1
            while True:
                my_row = BedRow.BedRow(line.strip())
                my_id = f"{my_row.id_str}.p{x}"
                row_len = my_row.get_exon_total_length()
                if cds_dict.get(my_id):
                    # Exclude reverse-strand predictions
                    # Note that in the orfipy output +/- is relative to the strand of the input transcript
                    # -> + always corresponds the strand the transcript was originally aligned to
                    if cds_dict.get(my_id)["strand"] != "+":
                        x += 1
                        continue
                    try:
                        if cds_dict[my_id]["start"] > 1:
                            # my_row.trim_upstream(cds_dict[my_id]['start'])
                            my_row.trim_bp_upstream(cds_dict[my_id]["start"])
                        # my_row.trim_downstream(cds_dict[my_id]['stop'] - cds_dict[my_id]['start'])
                        my_row.trim_bp_downstream(row_len - cds_dict[my_id]["stop"])

                        # Orfipy-style ID string
                        my_row.id_str = f"{my_row.id_str}.p{x}"
                        my_row.thick_start = my_row.start
                        my_row.thick_stop = my_row.thick_stop
                        export_rows.append(str(my_row) + "\n")
                        export_dict[my_row.id_str] = my_row
                    except Exception as e:
                        print(e)
                        dropped_count += 1
                else:
                    # No more predictions available -> break out of loop
                    break
                x += 1

    with open(f"{tmp_dir}/ExportCDS.bed", "w") as export_file:
        export_file.writelines(export_rows)

    # Now also include nested ORFs
    orf_for_blast_dict = dict()
    c = 0
    safe_drops = 0
    total = 0
    logger.info(f"Trimming {aligned_bed} to nested ORFs")
    with open(output_bed, "w") as export_orf_file:
        export_orf_file.writelines(export_rows)
        with open(nested_data, "r") as nested_file:
            total += 1
            for line in nested_file:
                cols = line.split("\t")
                # Extract ID
                # my_id = cols[0].split("_type")[0][1:]
                my_id = reformat_id_orfipy(cols[0][1:], False)
                orf = cols[0].split("_")[-1]
                offset = int(cols[1])

                # Entry is missing due to error in previous step
                if not export_dict.get(my_id):
                    if cols[0].find("(-)") != -1:
                        safe_drops += 1
                    else:
                        c += 1
                    # print(f"Missing: {my_id}")
                    continue

                # Deep copy in order to not alter the original record
                my_row = copy.deepcopy(export_dict[my_id])

                # The original, aligned ORF
                if offset == 0:
                    # Already written?
                    if not orf_for_blast_dict.get(my_row.id_str):
                        if my_row.stop < my_row.start:
                            print(f"Malformed: {my_id}")
                            c += 1
                        else:
                            orf_for_blast_dict[my_row.id_str] = my_row
                            export_orf_file.write(f"{str(my_row)}\n")
                    continue

                # Now modify according to ORF
                my_row.id_str = f"{my_row.id_str}_{orf}"

                # Offset is in bp now
                if offset * 3 + 1 < my_row.get_exon_total_length():
                    my_row.trim_bp_upstream(offset * 3 + 1)
                else:
                    continue

                my_row.stop += 1
                my_row.thick_stop = my_row.stop

                if my_row.strand == "-":
                    my_row.block_lens[-1] += 1
                else:
                    my_row.start -= 1
                    my_row.block_starts = [my_row.block_starts[0]] + [
                        x + 1 for x in my_row.block_starts[1:]
                    ]
                    my_row.stop -= 1

                if (
                    my_row.start + my_row.block_starts[-1] + my_row.block_lens[-1]
                    > my_row.stop
                ):
                    c += 1
                    # print(f"Blocks too long {my_id} / {my_row.block_starts[-1]}")
                    continue

                my_row.thick_start = my_row.start
                orf_for_blast_dict[my_row.id_str] = my_row
                if my_row.stop > my_row.start:
                    export_orf_file.write(f"{str(my_row)}\n")
                else:
                    c += 1

    if c > 0:
        logger.warning(f"Discarded {c} rows (safely)")


def export_cds_orfs(gff_file, aligned_bed, nested_data, output_bed, tmp_dir):
    """
    DEPRECATED
    Export CDS and nested ORFs from a GFF-based output of an ORF finder -> specifically Transdecoder
    :param gff_file: GFF ORF file -> chroms are transcript IDs
    :param aligned_bed: The transcript alignmend BED file
    :param nested_data: .tsv nested ORF file
    :param output_bed: Output path
    :param tmp_dir: Unified directory for temporary files to use
    :return:
    """

    cds_dict = dict()
    curr_end = 0

    # Load CDS predicted by Transdecoder
    logger.info(f"Loading Transdecoder CDS from {gff_file}")
    with open(gff_file, "r") as gff_file:
        for line in gff_file:
            if len(line) < 2:
                continue
            cols = line.split("\t")
            if cols[2] == "gene":
                curr_end = int(cols[4])
            if cols[2] == "CDS":
                my_id = cols[8].split(";")[0][7:]
                cds_dict[my_id] = {
                    "start": int(cols[3]),
                    "stop": int(cols[4]),
                    "final_stop": curr_end,
                    "strand": cols[6],
                }

    export_rows = []
    export_dict = dict()
    dropped_count = 0
    logger.info(f"Trimming {aligned_bed} to CDS coordinates")
    with open(aligned_bed, "r") as ref_file:
        for line in ref_file:
            x = 1
            while True:
                my_row = BedRow.BedRow(line.strip())
                my_id = f"{my_row.id_str}.p{x}"
                if cds_dict.get(my_id):
                    if cds_dict.get(my_id)["strand"] != "+":
                        x += 1
                        continue
                    try:
                        if cds_dict.get(my_id)["strand"] == "+":
                            if cds_dict[my_id]["start"] > 1:
                                my_row.trim_upstream(cds_dict[my_id]["start"] - 1)
                            my_row.trim_downstream(
                                cds_dict[my_id]["stop"] - cds_dict[my_id]["start"] + 1
                            )
                        my_row.id_str = f"{my_row.id_str}.p{x}"
                        # my_row.thick_start = my_row.start
                        # my_row.thick_stop = my_row.stop
                        export_rows.append(str(my_row) + "\n")
                        export_dict[my_row.id_str] = my_row
                    except Exception as e:
                        dropped_count += 1
                else:
                    # print(f"Not found {my_id}")
                    break
                x += 1

    with open(f"{tmp_dir}/ExportCDS.bed", "w") as export_file:
        export_file.writelines(export_rows)

    orf_for_blast_dict = dict()
    c = 0
    logger.info(f"Trimming {aligned_bed} to nested ORFs")
    with open(output_bed, "w") as export_orf_file:
        export_orf_file.writelines(export_rows)
        with open(nested_data, "r") as nested_file:
            for line in nested_file:
                cols = line.split("\t")
                # Extract ID
                my_id = cols[0].split("_type")[0][1:]
                orf = cols[0].split("_")[-1]
                offset = int(cols[1])

                # Entry is missing due to error in previous step
                if not export_dict.get(my_id):
                    continue

                # Deep copy in order to not alter the original record
                my_row = copy.deepcopy(export_dict[my_id])

                # The original, aligned ORF
                if offset == 0:
                    if not orf_for_blast_dict.get(my_row.id_str):
                        if my_row.stop < my_row.start:
                            c += 1
                        else:
                            orf_for_blast_dict[my_row.id_str] = my_row
                            export_orf_file.write(f"{str(my_row)}\n")
                    continue

                # Now modify according to ORF
                my_row.id_str = f"{my_row.id_str}_{orf}"

                # Offset is in bp now
                if offset * 3 + 1 < my_row.get_exon_total_length():
                    my_row.trim_upstream(offset * 3 + 1)
                else:
                    continue

                my_row.stop += 1
                my_row.thick_stop = my_row.stop

                if my_row.strand == "-":
                    my_row.block_lens[-1] += 1
                else:
                    my_row.start -= 1
                    my_row.block_starts = [my_row.block_starts[0]] + [
                        x + 1 for x in my_row.block_starts[1:]
                    ]
                    my_row.stop -= 1

                if (
                    my_row.start + my_row.block_starts[-1] + my_row.block_lens[-1]
                    > my_row.stop
                ):
                    c += 1
                    continue

                my_row.thick_start = my_row.start
                orf_for_blast_dict[my_row.id_str] = my_row
                if my_row.stop > my_row.start:
                    export_orf_file.write(f"{str(my_row)}\n")
                else:
                    c += 1
    logger.warning(f"Discarded {c} rows (safely)")


def get_nested_orfs(input_file):
    """
    Wrapper for splitAtMs script
    :param input_file: File path passed to splitAtMs
    :return:
    """
    process = subprocess.run(
        f"/projects/hillerlab/genome/bin/scripts/splitAtMs.perl -minLen 30 -minPercent 0.25 {input_file}",
        shell=True,
        check=True,
    )
    return process.returncode == 0


def blast_main(
    tmp_dir,
    transcript_fasta_file,
    transcript_aligned_bed,
    blast_db,
    use_blast=False,
    use_nextflow=False,
    use_orfipy=True,
):
    # Step 1: Run Transdecoder on input data
    use_diamond = not use_blast

    if not use_orfipy:
        logger.info("Running Transdecoder")
        try:
            run_transdecoder_long_orfs(transcript_fasta_file, f"{tmp_dir}/ORF")
        except Exception as e:
            logger.error(f"Transdecoder failed {e}")
            exit(-1)

        logger.info("Transdecoder finished successfully")
        if not os.path.exists(f"{tmp_dir}/ORF/longest_orfs.pep"):
            logger.error(
                f"Failed to find Transdecoder results at: {tmp_dir}/ORF/longest_orfs.pep"
            )
            exit(-1)

        logger.info("Converting Transdecoder results to BED for later")
        subprocess.run(
            f"gff2bed < {tmp_dir}/ORF/longest_orfs.gff3 > {tmp_dir}/ORF/longest_orfs.bed",
            shell=True,
        )
    else:
        try:
            logger.info("Running orfipy")
            run_orfipy(transcript_fasta_file, f"{tmp_dir}/ORF")
        except Exception as e:
            logger.error(f"Orfipy failed {e}")
            exit(-1)

        logger.info("Orfipy finished successfully")

    # Step 2: Get nested ORFs -> Transdecoder/longest_orfs.pep.nestedORFs.fa
    if not use_orfipy:
        get_nested_orfs(f"{tmp_dir}/ORF/longest_orfs.pep")
        logger.info(
            f"Sequences to BLAST written to {tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa"
        )
    else:
        subprocess.run(
            f"cat {tmp_dir}/ORF/orfs_raw.pep | linearize_fasta > {tmp_dir}/ORF/orfs.pep",
            shell=True,
        )
        get_nested_orfs(f"{tmp_dir}/ORF/orfs.pep")
        shutil.move(
            f"{tmp_dir}/ORF/orfs.pep.nestedORFs.fa",
            f"{tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa",
        )
        logger.info(
            f"Sequences to BLAST written to {tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa"
        )

    # TODO: ORF + length cutoffs

    logger.info(
        "Reducing number of sequences to BLAST by multiplexing common sequences"
    )
    if not use_orfipy:
        multiplex_dataset(
            f"{tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa",
            f"{tmp_dir}/ORF/multiplexed.fa",
            f"{tmp_dir}/ORF/multiplexed.index",
        )
    else:
        multiplex_dataset(
            f"{tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa",
            f"{tmp_dir}/ORF/multiplexed.fa",
            f"{tmp_dir}/ORF/multiplexed.index",
        )

    # Step 3: BLAST P
    logger.info("Running BLAST")
    if use_diamond:
        result = run_diamond(
            f"{tmp_dir}/ORF/multiplexed.fa",
            f"{tmp_dir}/BLASTP/multiplexed.fa.fmt6",
            blast_db,
            tmp_dir,
            chunk_size=3000,
            nextflow=use_nextflow,
        )
    else:
        result = run_blastp(
            f"{tmp_dir}/ORF/multiplexed.fa",
            f"{tmp_dir}/BLASTP/multiplexed.fa.fmt6",
            blast_db,
            tmp_dir,
            chunk_size=500,
        )

    if result.returncode != 0:
        logger.error("Error running BLAST, dumping output:")
        print(result.stdout, "========", result.stderr)
        exit(-1)

    # Inflate and sort BLAST output
    de_multiplex_blast(
        f"{tmp_dir}/ORF/multiplexed.index",
        f"{tmp_dir}/BLASTP/multiplexed.fa.fmt6",
        f"{tmp_dir}/BLASTP/longest_orfs.pep.nestedORFs.fa.unsorted.fmt6",
    )
    subprocess.run(
        f"sort -k1 < {tmp_dir}/BLASTP/longest_orfs.pep.nestedORFs.fa.unsorted.fmt6 > {tmp_dir}/BLASTP/longest_orfs.pep.nestedORFs.fa.fmt6",
        shell=True,
    )

    # Step 4: producing rich tracks
    logger.info("Creating rich tracks for ORFs")
    if not use_orfipy:
        export_cds_orfs(
            f"{tmp_dir}/ORF/longest_orfs.gff3",
            transcript_aligned_bed,
            f"{tmp_dir}/ORF/longest_orfs.pep.nestedORFs.tsv",
            f"{tmp_dir}/BLASTP/trimmed_alignments.bed",
            tmp_dir,
        )
    else:
        export_cds_orfs_bed(
            f"{tmp_dir}/ORF/orfs.bed",
            transcript_aligned_bed,
            f"{tmp_dir}/ORF/orfs.pep.nestedORFs.tsv",
            f"{tmp_dir}/BLASTP/trimmed_alignments.bed",
            tmp_dir,
        )
