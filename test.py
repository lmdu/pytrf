import pytrf

s = "AAGAAGAAGAAGCCGAGAAGGTAGATAG"

atrs = pytrf.ATRFinder('s1',s).as_list()

print(atrs)

