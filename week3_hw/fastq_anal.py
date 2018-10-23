import sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import numpy as np
import seaborn as sns
from random import choice

filename = "mytestfile.fastq"

fastq_file = open(filename)
lines = fastq_file.read().strip().split("\n")
fastq_file.close()


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

df = pd.DataFrame(data=seq_probs)

summary_df = pd.DataFrame(data={"means":means,"medians":medians,"mins":mins,"maxs":maxs})


f2, axs = plt.subplots(2, 1, figsize=(10, 8))
f2.suptitle("Histogram of Min & Max probabilities of Nucleotide Correctness over 100 sequences",fontsize=15)
axs[0].hist(summary_df["mins"],bins=100)
axs[0].set_xlabel("Minimum of Probabilities")
axs[1].hist(summary_df["maxs"],range=[0.999,1.0],bins=100)
axs[1].set_xlabel("Maximum of Probabilities")

plt.show()

range_means = max(summary_df["means"])-min(summary_df["means"])
range_medians = max(summary_df["medians"])-min(summary_df["medians"])

f1, ax = plt.subplots(1, 1, figsize=(10, 8))

ax.set_xlim(min(summary_df["means"])-range_means/4, max(summary_df["means"])+range_means/4)
ax.set_ylim(min(summary_df["medians"])-range_medians/4, max(summary_df["medians"])+range_medians/4)
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)

current_palette = sns.color_palette("GnBu_d", n_colors=5)
cmap = ListedColormap(sns.color_palette(current_palette).as_hex())

def onpick3(event):
    """creates a violin plot based on the individual nucleotide accuracies for the index at the clicked point"""
    ind = event.ind
    print('onpick3 scatter:', ind, np.take(means, ind), np.take(medians, ind))
    i = choice(ind)
    f3,ax3 = plt.subplots(1, 1, figsize=(5, 5))
    ax3.violinplot(seq_probs[i], vert=False)
    ax3.set_title("Violin Plot of Nucleotide Accuracy for Seq. {}".format(i))
    ax3.set_xlabel("Accuracy of Nucleotide")
    plt.show()

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--',c="red")

    

ax.scatter(summary_df["means"],summary_df["medians"],c=range(len(means)),cmap=cmap,picker=True)
f1.canvas.mpl_connect('pick_event', onpick3)
# if the probabilities follow a symmetric distribution the mean & median should be pretty close to each other
# the scatter plot should follow the line y=x
abline(1,0)
plt.xlabel("Avg Prob of Nucleotide Correctness",fontsize=12)
plt.ylabel("Median Prob of Nucleotide Correctness",fontsize=12)
plt.title("Mean vs. Median of Probability of Nucleotide Correctness over 100 sequences",fontsize=15)
plt.show()




