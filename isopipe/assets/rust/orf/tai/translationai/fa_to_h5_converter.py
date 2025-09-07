########################################################################################
"""This program takes as input the .fa file and outputs a .h5 file with datapoints
of the form (X, Y), which can be understood by Keras models."""
# Usage: fa_to_h5_converter.py input.fa output.h5
########################################################################################

import h5py
import sys
import time
import re
import numpy as np
from translationai.utils import *


def parse_fasta_chunks(fnIn, pad_len, chunk_size):
    """
    Generator that yields (chunk_size)-sized chunks of padded sequences and metadata.
    Each yielded chunk is a list of tuples:
    (padded_seq, strand, tx_start, tx_end, TIS_list, TTS_list)
    """
    chunk = []
    with open(fnIn, "r") as file:
        header_data = None
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                words = line[1:].split(":")
                strand = "+"
                try:
                    words2 = re.split(r"\(|\)", words[1])
                except Exception:
                    raise Exception(f"Seq header format error: {line}")
                if len(words2) not in {5, 7}:
                    raise Exception(f".fa header format error!\n{line}")
                try:
                    TIS_pos, TTS_pos = words2[-2].split(" ")
                except ValueError:
                    raise Exception(f".fa TIS/TTS format error!\n{line}")
                header_data = (
                    strand,
                    "0",
                    None,
                    [TIS_pos],
                    [TTS_pos],
                )  # TX_END is filled later
            else:
                seq = line.upper()
                if not re.match(r"^[ACGTN]+$", seq):
                    raise ValueError(
                        f"Invalid sequence: {seq}\nOnly A, C, G, T, and N are allowed."
                    )
                padded_seq = "N" * pad_len + seq + "N" * pad_len
                seq_len = len(seq)
                strand, tx_start, _, tis, tts = header_data
                tx_end = str(seq_len - 1)
                chunk.append((padded_seq, strand, tx_start, tx_end, tis, tts))
                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
        if chunk:
            yield chunk


def main(fnIn, fnOut, CHUNK_SIZE=1):
    start_time = time.time()

    print("----Reading and converting .fa file to .h5 format for translationAI----")
    print(f"    *Input: {fnIn}")
    print(f"    *Output: {fnOut}")
    print(f"    *Chunk size: {CHUNK_SIZE}")

    pad_len = CL_max // 2  # from translationai.utils

    h5f2 = h5py.File(fnOut, "w")

    for i, chunk in enumerate(parse_fasta_chunks(fnIn, pad_len, CHUNK_SIZE)):
        X_batch = []
        Y_batch = [[] for _ in range(1)]  # Supports multiple Y outputs

        for padded_seq, strand, tx_start, tx_end, tis, tts in chunk:
            X, Y = create_datapoints(padded_seq, strand, tx_start, tx_end, tis, tts)
            X_batch.extend(X)
            for t in range(1):
                Y_batch[t].extend(Y[t])

        X_batch = np.asarray(X_batch, dtype="int8")
        for t in range(1):
            Y_batch[t] = np.asarray(Y_batch[t], dtype="int8")

        h5f2.create_dataset("X" + str(i), data=X_batch)
        h5f2.create_dataset("Y" + str(i), data=Y_batch)

        # print(f'    *Wrote block {i} ({len(chunk)} sequences, {len(X_batch)} datapoints)')

    h5f2.close()
    elapsed = int(time.time() - start_time)
    print(f"----Conversion completed successfully! Time used: {elapsed} seconds.----")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fa_to_h5_converter.py input.fa output.h5")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
