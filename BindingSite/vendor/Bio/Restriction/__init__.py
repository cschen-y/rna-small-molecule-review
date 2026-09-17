








"""Restriction Digest Enzymes.

Examples
--------

.. code-block:: pycon

    >>> from Bio.Seq import Seq
    >>> from Bio.Restriction import *
    >>> pBs_mcs = "GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTG"
    >>> pBs_mcs += "CAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGCCACCGCGGTGGAGCTC"
    >>> seq = Seq(pBs_mcs)  # Multiple-cloning site of pBluescript SK(-)
    >>> a = Analysis(AllEnzymes, seq)
    >>> a.print_that()  # no argument -> print all the results
    AbaSI      :  10, 12, 13, 16, 17, 18, 19, 20, 22, 23, 24, 25, 25, 26, 27...
    BmeDI      :  7, 8, 8, 9, 9, 13, 14, 15, 16, 17, 18, 19, 19, 21, 21...
    YkrI       :  10, 12, 13, 16, 16, 17, 19, 20, 21, 22, 23, 24, 25, 25, 26...

    BmeDI      :  1, 2, 7, 8, 8, 9, 9, 13, 14, 15, 16, 17, 18, 19...
    AccII      :  98.
    AciI       :  86, 90, 96, 98...

    Enzymes which do not cut the sequence.

    AspLEI    BstHHI    CfoI      CviAII    FaeI      FaiI      FatI      GlaI
    HhaI      Hin1II    Hin6I     HinP1I    HpyCH4IV  HpySE526I Hsp92II   HspAI
    MaeII     MseI      NlaIII    SaqAI     TaiI      Tru1I     Tru9I...
    <BLANKLINE>
    >>> b = a.blunt()  # Analysis with blunt enzymes
    >>> a.print_that(b)  # Print results for blunt cutters
    AccII      :  98.
    AfaI       :  4.
    AluBI      :  40, 106.
    AluI       :  40, 106.
    Bsh1236I   :  98.
    BshFI      :  10, 89.
    BsnI       :  10, 89.
    BspANI     :  10, 89...

    Enzymes which do not cut the sequence.

    FaiI      GlaI      CdiI      MlyI      SchI      SspD5I    AanI...
    <BLANKLINE>

"""  

from Bio.Restriction.Restriction import *  








































































































































