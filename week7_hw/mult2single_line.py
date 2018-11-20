import sys
import re

unknown_gap = r'N{10,}'

cur = ""
contigs = []
for line in sys.stdin:
    line = line.strip()
    if line.startswith(">"):
        if len(cur)>0:
            contigs.extend(re.split(unknown_gap,cur))
            cur = ""
    else:
        cur+=line
if len(cur)>0:
    contigs.extend(re.split(unknown_gap,cur))
    cur = ""


out = "\n".join(contigs)

print(out)