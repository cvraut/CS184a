import sys

cur = ""
contig_lens = []
for line in sys.stdin:
    line = line.strip()
    contig_lens.append(len(line))

out = " ".join(map(str,contig_lens))

print(out)