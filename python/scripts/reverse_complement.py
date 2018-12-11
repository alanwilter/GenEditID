from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

seq = Seq('TGGACAGCAACAAGGGTCAA', IUPAC.unambiguous_dna)
seq = str(seq.reverse_complement())
print(seq)
