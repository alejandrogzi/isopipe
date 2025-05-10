from typing_extensions import Optional, Tuple


class ColumnCountError(Exception):
    """
    Raised when the column count of a BED entry does not correspond to a valid format.
    """

    def __init__(self, count):
        super().__init__(
            f"Found {count} tab-separated columns, expected one of [6, 8, 9, 12]"
        )
        pass


class MismatchedBlockCounts(Exception):
    """
    Raised when the block indices and block sizes are of different length.
    """

    def __init__(self, indices_length, block_size_length, block_count):
        super().__init__(
            f"The total number of block indices ({indices_length}) is inconsisten to the total number of block sizes ({block_size_length}) or the block count {block_count}."
        )
        pass


class OutOfExonError(Exception):
    """ """

    def __init__(self, genomic_position, row_id):
        super().__init__(
            f"The genomic position ({genomic_position}) is not inside a block of {row_id}"
        )
        pass


class BedRow:
    """
    A Class representing a single row in a BED file.
    """

    def __init__(self, row: str):
        features = row.strip().split("\t")
        self.bed_x = len(features)

        if self.bed_x not in [6, 8, 9, 12]:
            raise ColumnCountError(self.bed_x)

        # BED6
        self.chrom = features[0]
        self.start = int(features[1])
        self.stop = int(features[2])
        self.id_str = features[3]
        self.score = features[4]
        self.strand = features[5]

        if self.bed_x > 6:
            self.thick_start = int(features[6])
            self.thick_stop = int(features[7])
        if self.bed_x > 7:
            self.rgb = features[8]
        if self.bed_x > 8:
            self.block_count = int(features[9])

            # There are some inconsistencies in the usage of a trailing comma in BED files depending on source
            if features[10].endswith(","):
                self.block_lens = [int(x) for x in features[10].split(",")[:-1]]
            else:
                self.block_lens = [int(x) for x in features[10].split(",")]
            # It may be inconsistent across block sizes / starts
            if features[11].endswith(","):
                self.block_starts = [int(x) for x in features[11].split(",")[:-1]]
            else:
                self.block_starts = [int(x) for x in features[11].split(",")]

    def __str__(self):
        if self.bed_x == 6:
            return "\t".join(
                str(x)
                for x in [
                    self.chrom,
                    str(self.start),
                    str(self.stop),
                    self.id_str,
                    self.score,
                    self.strand,
                ]
            )
        elif self.bed_x == 8:
            return "\t".join(
                str(x)
                for x in [
                    self.chrom,
                    str(self.start),
                    str(self.stop),
                    self.id_str,
                    self.score,
                    self.strand,
                    self.thick_start,
                    self.thick_stop,
                ]
            )
        elif self.bed_x == 9:
            return "\t".join(
                str(x)
                for x in [
                    self.chrom,
                    str(self.start),
                    str(self.stop),
                    self.id_str,
                    self.score,
                    self.strand,
                    self.thick_start,
                    self.thick_stop,
                    self.rgb,
                ]
            )
        elif self.bed_x == 12:
            return "\t".join(
                str(x)
                for x in [
                    self.chrom,
                    str(self.start),
                    self.stop,
                    self.id_str,
                    self.score,
                    self.strand,
                    self.thick_start,
                    self.thick_stop,
                    self.rgb,
                    self.block_count,
                    str(self.block_lens)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(" ", ""),
                    str(self.block_starts)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(" ", ""),
                ]
            )

    def get_exon_total_length(self):
        """
        Returns the sum of all block
        :return:
        """
        return sum(self.block_lens)

    def get_total_block_length(self):
        """
        Returns the sum of all block
        :return:
        """
        return sum(self.block_lens)

    def is_genomic_position_in_block(
        self, genomic_position: int, force_exon=False
    ) -> bool:
        """
        Safely wrap genomic position <-> block position conversion
        :param force_exon: If set, enforce that the block position must also be within the thickStart / Stop bounds.
        :param genomic_position:
        :return: True if the position is within a block
        """
        if force_exon and (
            genomic_position < self.thick_start or genomic_position > self.thick_stop
        ):
            return False

        if self.start > genomic_position or self.stop < genomic_position:
            return False

        # In BED6 everything between start/stop is in a block
        if self.bed_x == 6:
            return True

        try:
            self.genomic_position_to_block_position(genomic_position)
            return True
        except OutOfExonError:
            return False

    def get_bp_pos_block(self, pos: int) -> Tuple[int, int]:
        """
        For a given base pair position, find the block it is located in, and it's position in the block. The bp position
        is always in (+) direction, as is in the BED standard.
        :param pos:
        :return:
        """

        curr_pos = 0
        for x, block in enumerate(self.block_lens):
            # Advance to the next block
            if curr_pos + block < pos:
                curr_pos += block
            # The position is in the current block
            else:
                # BP left in block
                remainder = pos - curr_pos
                return x, remainder

        return -1, -1

    def genomic_position_to_block_position(
        self, genomic_pos: int
    ) -> Optional[Tuple[int, int]]:
        """
        Convert a genomic position into the internal reference bp count of the BedRow.
        :param genomic_pos:
        :return:
        """
        # Out of feature
        if genomic_pos < self.start or genomic_pos > self.stop:
            raise OutOfExonError(genomic_pos, self.id_str)

        for x, block in enumerate(self.block_lens):
            # The end of the current block is downstream of the genomic position
            if self.start + self.block_starts[x] + block < genomic_pos:
                continue
            else:
                delta = genomic_pos - self.start - self.block_starts[x]
                if delta < 0:
                    raise OutOfExonError(genomic_pos, self.id_str)
                else:
                    return x, delta

    def block_position_to_genomic_position(self, bp_pos: int) -> int:
        """
        Convert a block position in (+) direction into a genomic position.
        :param bp_pos:
        :return:
        """
        block, remainder = self.get_bp_pos_block(bp_pos)
        return self.start + self.block_starts[block] + remainder

    def smooth_zero_leading_block(self):
        """
        Invoke this when the leading block has a length of 0
        :return:
        """

        # Remove the starting block
        self.block_lens = self.block_lens[1:]
        self.block_starts = self.block_starts[1:]
        leading_block = self.block_starts[0]

        # Move all blocks
        for x, block in enumerate(self.block_starts):
            self.block_starts[x] = block - leading_block

        # Move the start
        self.start = self.start + leading_block
        # And thick starts if necessary
        if self.thick_start < self.start:
            self.thick_start = self.start

    def trim_bp_downstream(self, bp: int, allow_invert=True):
        """
        Trim a set amount of base pairs off the downstream end. This will remove blocks and adjust thick borders if necessary.
        :param allow_invert: Suppress automatic inversion to upstream trim if on (-) strand
        :param bp:
        :return:
        """
        if bp <= 0:
            return

        if self.strand == "-" and allow_invert:
            return self.trim_bp_upstream(bp, False)

        # Find BP position
        block_index, remainder = self.get_bp_pos_block(
            self.get_exon_total_length() - bp
        )

        # Remove blocks if necessary
        self.block_starts = self.block_starts[: block_index + 1]
        self.block_lens = self.block_lens[: block_index + 1]
        self.block_count = block_index + 1

        # Trim the last block
        if remainder != self.block_lens[-1]:
            self.block_lens[-1] = remainder

        # Adjust stop
        self.stop = self.start + self.block_starts[-1] + self.block_lens[-1]
        if self.thick_stop > self.stop:
            self.thick_stop = self.stop

    def trim_bp_upstream(self, bp: int, allow_invert=True):
        """
        Trim a set amount of base pairs off the upstream end. This will remove blocks and adjust thick borders if necessary.
        :param allow_invert: Suppress automatic inversion to downstream trim if on (-) strand
        :param bp:
        :return:
        """
        if bp <= 0:
            return

        if self.strand == "-" and allow_invert:
            return self.trim_bp_downstream(bp, False)

        # Find BP position
        block_index, remainder = self.get_bp_pos_block(bp)

        # Remove blocks if necessary
        self.block_starts = self.block_starts[block_index:]
        self.block_lens = self.block_lens[block_index:]
        self.block_count = self.block_count - block_index

        # Shorten the block, move the start
        self.block_lens[0] = self.block_lens[0] - remainder
        self.block_starts[0] = self.block_starts[0] + remainder

        # Adjust start
        shift = self.block_starts[0]
        self.start = self.start + self.block_starts[0]
        if self.thick_start < self.start:
            self.thick_start = self.start

        self.block_starts = [x - shift for x in self.block_starts]

        if self.block_lens[0] == 0:
            self.smooth_zero_leading_block()

    def trim_to_genomic_position_downstream(self, genomic_position):
        """
        Takes a genomic position and trims everything downstream (in strand direction) of that position.
        :param genomic_position:
        :return:
        """

        # Find BP position
        cut = self.genomic_position_to_block_position(genomic_position)
        # Convert into trimmable block position
        walked = sum(self.block_lens[: cut[0]]) + cut[1]
        # Trim
        if self.strand == "+":
            self.trim_bp_downstream(self.get_exon_total_length() - walked)
        else:
            self.trim_bp_downstream(walked)

    def trim_to_genomic_position_upstream(self, genomic_position):
        """
        Takes a genomic position and trims everything upstream (in strand direction) of that position.
        :param genomic_position:
        :return:
        """
        # Find BP position
        cut = self.genomic_position_to_block_position(genomic_position)
        # Convert into trimmable block position
        walked = sum(self.block_lens[: cut[0]]) + cut[1]
        # Trim
        if self.strand == "-":
            self.trim_bp_upstream(self.get_exon_total_length() - walked)
        elif self.strand == "+":
            self.trim_bp_upstream(walked)

    def genomic_position_to_sequence_position(self, genomic_position: int) -> int:
        """
        Convert a given genomic position to a sequence position
        :param genomic_position:
        :return:
        """
        # Internal reference
        cut = self.genomic_position_to_block_position(genomic_position)
        walked = sum(self.block_lens[: cut[0]]) + cut[1]

        # The walked bp is always in (+) -> invert if on (-) strand
        if self.strand == "-":
            return self.get_total_block_length() - walked
        else:
            return walked

    def sequence_position_to_genomic_position(self, sequence_position: int) -> int:
        """
        Convert a sequence position (with the sequence spanning all blocks) to a genomic position
        :param sequence_position:
        :return:
        """

        # First, convert the sequence position into a block position
        if self.strand == "+":
            block_pos = sequence_position
        else:
            block_pos = self.get_total_block_length() - sequence_position

        return self.block_position_to_genomic_position(block_pos)

    ## Deprecated methods

    def trim_upstream(self, pos: int, thick=False, smoothing=False):
        raise DeprecationWarning("This function was deprecated in March 2025")

    def trim_downstream(self, pos: int, thick=False):
        raise DeprecationWarning("This function was deprecated in March 2025")
