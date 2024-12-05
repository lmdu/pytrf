import pytrf

#seq = "ATTGAATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCATCAAATGGAATCGAATGGAATCAACATCAAATGGAATCTAATGGAATCATTGAACGGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCAAATGGAATGGAATGGAAT"
#seq = "TCCTCCTCCTCCTCTTTCCATTAAGTT"
#seq = "TGTGTGTGTGTGTGTTTTGTTTTTCAACTTAAAAA"
#seq = "AAAAAAAAAAAAAAAAGAAATCCAAATAAAATTT"
seq = "GTTTTGTTTTGTTTTGTTTTTTGAGACAGAGTTTC"

for atr in pytrf.ATRFinder('chr1', seq):
	print(atr.as_string())
