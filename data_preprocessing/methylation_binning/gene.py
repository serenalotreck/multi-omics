"""
Gene class definition to use for binning methylation.

Author: Serena G. Lotreck
"""

class Gene:

    def __init__(self, gff_line):
        """
        Instantiate an instance of the Gene class. 

        parameters:
            gff_line, list of str: line from a gff file split on \t that
                represents the gene

        returns: None
        """
        self.seqid = gff_line[0]
        self.start = int(gff_line[3])
        self.end = int(gff_line[4])

        self.CG = []
        self.CHG = []
        self.CHH = []

    def 
