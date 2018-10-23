import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filename = "mytestfile.fastq"

fastq_file = open(filename)

lines = fastq_file.read().strip().split("\n")

seq_probs = []
means = []
medians = []
mins = []
maxs = []
for i,line in enumerate(lines):
    if i%4 == 3:
        seq_probs.append([1-10.0**(-1*(ord(c)-33)/10.0) for c in line.strip()])
        means.append(np.mean(seq_probs[-1]))
        medians.append(np.median(seq_probs[-1]))
        mins.append(min(seq_probs[-1]))
        maxs.append(max(seq_probs[-1]))
df = pd.DataFrame(seq_probs)

plt.scatter(means,medians)
plt.show()

fastq_file.close()

