#!/usr/bin/env python3
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.utils import split_nth
from singlecellmultiomics.tagtools import tagtools


def fromTideHunterBam(self, pysamRecord):
    # TideHunter v1.4.1
    # Each samfile reads contains 1 TideHunterFragment
    # All of the following tag except CL are unique to each TideHunter Fragment.
    # TF: TideHunterFragment ID
    # CL: Cyclomic full length (CL is shared with all fragments with same readname-- same nanopre read.)
    # TL: Tandem repeat length (tandem repeat length including backbone or not base on TB tag)
    # TB: If 1: forward complete(supplied part) backbone found around the fragment, 2: reverse found, 0: not found
    # CS: Consensus score for the TF(average matching score for each repeat matched to the consensus 0-100%)
    # TC: Tandem repeat count. How many tandem repeats are found in this molecule. (Similar to NH tag)
    # TS: SubPositions of the TideHunter Fragments repeats (transfer to 0-based). [start, first, last]
    ordered_keys = ['readname',
                    "TF",
                    "CL",
                    "start",
                    "end",
                    "TL",
                    "TC",
                    "CS",
                    "TB",
                    "first of subPos"]
    # start enumerating from 1
    ordered_values = pysamRecord.query_name.strip().split('_')
    keyValues = zip(ordered_keys, ordered_values)
    ts = []
    # skip read name
    self.addTagByTag("LY", "")
    for (key, value) in keyValues:
        if len(key) == 2:
            self.addTagByTag(key, value, isPhred=False)
        else:
            if key == 'readname':
                pass
            else:
                ts.append(value)
                self.addTagByTag("TS", tuple(ts), isPhred=False)


class TideHunterFlagger(DigestFlagger):
    """This query name flagger converts values between "underline" to tags"""

    def __init__(self, block_assignments, **kwargs):
        """Initialise TideHunterFlagger

        Args:
            block_assignments(list) : list of two letter codes to assign blocks to

        """
        self.block_assignments = block_assignments
        self.origin_fields =   # amount of ':' in original read name
        # Verify if all of the assignments are 2 letters:
        if not all((len(b) == 2 for b in block_assignments)):
            for b in block_assignments:
                if len(b) != 2:
                    raise ValueError(f'Tag {b} is not two letters long')

        DigestFlagger.__init__(self, **kwargs)

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue

            # Split original read name from added data
            origin, rest = split_nth(read.query_name, ':', self.origin_colons)
            # Reset the read name
            read.query_name = origin
            # Write the tags
            for tag, value in zip(self.block_assignments, rest.split(':')):
                read.set_tag(tag, value)


