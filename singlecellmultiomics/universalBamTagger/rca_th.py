#!/usr/bin/env python3
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.utils import split_nth


class RCA_Tidehunter_Flagger(DigestFlagger):
    def __init__(self, separator=None, origin_colons=None, SMTag = None, **kwargs):
        self.separator = "_"
        self.SMTag = SMTag
        # separator count in filename. 0: there should be no "_" in the filename
        self.origin_colons = 1
        DigestFlagger.__init__(self, **kwargs)

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue
            # Split original read name from added data, return 2 strings
            origin, rest = split_nth(read.query_name, self.separator, self.origin_colons)
            # set query name to the original nanopore read name
            read.query_name = origin
            # add original read name as a tag
            read.set_tag("SM", self.SMTag)
            read.set_tag('CR', origin)
            # select useful added data and add tags
            repN, read_fl, start, end, tl, tc_copies, c_score, dir_td, sub_po  = rest.split(self.separator)
            read.set_tag('TF', repN)
            read.set_tag('CL', read_fl)
            read.set_tag('TL', tl)
            read.set_tag('TD', int(dir_td))
            read.set_tag('TC', float(tc_copies))
            read.set_tag('CS', float(c_score))
            read.set_tag('TS', [int(start), int(sub_po.split(',')[0]),int(end)])
