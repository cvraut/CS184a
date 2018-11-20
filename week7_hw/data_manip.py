import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from collections import Counter
from bisect import bisect

ref_contig_cnts_fname = "cnts.ref.d_melanogaster.contig_cnts.out"
exp_contig_cnts_fname = "cnts.iso1_onp_a2_1kb_30x_LD0_K25_KCOV2_ADAPT0.01_MINOVL35_RMCHIM1_asm_new_np_pilon_pilon.contig_cnts.out"
#"cnts.iso1_onp_keep_bn_keep_ngs.contig_cnts.out"

n_bins = 15

def get_contig_list(cnts_fname):
    cnts_file = open(cnts_fname)
    cnts = map(int,cnts_file.read().strip().split())
    cnts_file.close()
    return list(filter(lambda x:x>0,cnts))

def filter_cnts(ctgs,thresh = 5e4):
    s = []
    b = []
    for i in range(len(ctgs)):
        if ctgs[i] < thresh:
            s.append(ctgs[i])
        else:
            b.append(ctgs[i])
    return s,b
    
ref_contig_cnts = get_contig_list(ref_contig_cnts_fname)
exp_contig_cnts = get_contig_list(exp_contig_cnts_fname)


f, ax = plt.subplots(1, 1, figsize=(10, 8))
f.suptitle("histogram of ref vs exp contig sizes",fontsize=15)
colors = ['red', 'lime']
leg_labels = ["ref","exp"]
ax.hist([ref_contig_cnts,exp_contig_cnts], n_bins, histtype='bar', color=colors, label=leg_labels)
ax.legend(prop={'size': 10})
ax.set_xlabel("contig lengths")
#plt.show()

s_ref_contig_cnts,b_ref_contig_cnts = filter_cnts(ref_contig_cnts)
s_exp_contig_cnts,b_exp_contig_cnts = filter_cnts(exp_contig_cnts)

smalls = [s_ref_contig_cnts,s_exp_contig_cnts]
print(len(smalls),len(smalls[0]),len(smalls[1]))
bigs = [b_ref_contig_cnts,b_exp_contig_cnts]
print(len(bigs),len(bigs[0]),len(bigs[1]))

s_hist, s_bins, _ = plt.hist(smalls, bins=n_bins)
s_bins[0]+=1e-10
s_logbins = np.logspace(np.log10(s_bins[0]),np.log10(s_bins[-1]),len(s_bins))
b_hist, b_bins, _ = plt.hist(bigs, bins = n_bins)
b_bins[0]+=1e-10
b_logbins = np.logspace(np.log10(b_bins[0]),np.log10(b_bins[-1]),len(b_bins))

plt.close()


f, axes = plt.subplots(2, 1, figsize=(10, 10),gridspec_kw = {'height_ratios':[4, 2]})
f.suptitle("Distribution of Reference and Experimental Contig sizes",fontsize=15)
axes[0].hist(smalls,bins=s_logbins,histtype='bar', color=colors, label=leg_labels)
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlabel("Contigs < 50KB")
axes[0].legend(prop={'size': 10})

axes[1].hist(bigs,bins=b_logbins,color=colors, label=leg_labels)
axes[1].set_xscale('log')
axes[1].set_xlabel("Contigs >= 50KB")
axes[1].legend(prop={'size': 10})
plt.show()


plt.close()

ref_tot = sum(ref_contig_cnts)
exp_tot = sum(exp_contig_cnts)



ref_contig_cnts = bigs[0]
exp_contig_cnts = bigs[1]

ref_contig_dict = Counter(ref_contig_cnts)
ref_contigs = sorted(ref_contig_dict.keys(),reverse=True)
cp_ref = []

for contig in ref_contigs:
    cur = 0
    if len(cp_ref)>0:
        cur+=cp_ref[-1]
    cp_ref.append((cur+contig*ref_contig_dict[contig]))

cp_ref= list(map(lambda x:x/ref_tot,cp_ref))


exp_contig_dict = Counter(exp_contig_cnts)
exp_contigs = sorted(exp_contig_dict.keys(),reverse=True)
cp_exp = []
for contig in exp_contigs:
    cur = 0
    if len(cp_exp)>0:
        cur+=cp_exp[-1]
    cp_exp.append((cur+contig*exp_contig_dict[contig]))

cp_exp= list(map(lambda x:x/exp_tot,cp_exp))
    
#ref_contig_cnts.sort(reverse=True)
#exp_contig_cnts.sort(reverse=True)

#cp_ref = np.cumsum(ref_contig_cnts)/sum(ref_contig_cnts)
#cp_exp = np.cumsum(exp_contig_cnts)/sum(exp_contig_cnts)

x_ref = np.linspace(min(ref_contigs),max(ref_contigs),100)
x_exp = np.linspace(min(exp_contigs),max(exp_contigs),100)

xs_ref = np.array(ref_contigs[::-1])
ys_ref = np.array(cp_ref[::-1])
spl_ref = UnivariateSpline(xs_ref,ys_ref,k=2)

xs_exp = np.array(exp_contigs[::-1])
ys_exp = np.array(cp_exp[::-1])
spl_exp = UnivariateSpline(xs_exp,ys_exp,k=2)

f, ax = plt.subplots(1, 1, figsize=(10, 8))

ax.scatter(ref_contigs,cp_ref,c=colors[0],alpha = 0.7,s=3)
ax.scatter(exp_contigs,cp_exp,c=colors[1],alpha = 0.7,s=3)

ax.plot(x_ref,spl_ref(x_ref),color = colors[0],lw=3,alpha=0.3,linestyle = "--",label=leg_labels[0]+" Smoothed")
ax.plot(x_exp,spl_exp(x_exp),color = colors[1],lw=3,alpha=0.3,linestyle = "--",label=leg_labels[1]+" Smoothed")

ax.plot(ref_contigs,cp_ref,c=colors[0],alpha = 0.85,label=leg_labels[0])
ax.plot(exp_contigs,cp_exp,c=colors[1],alpha = 0.85,label=leg_labels[1])




ax.set_xscale('log')
ax.invert_xaxis()

ax.legend(prop={'size': 10})


plt.xlabel("Contig size",fontsize=12)
plt.ylabel("Proportion of genome",fontsize=12)
plt.title("cdf of proportion of genome filled vs smallest contig added\nExcluding Contigs < 50KB",fontsize=15)
plt.show()
plt.close()

###################################################### smalls


ref_contig_cnts = smalls[0]
exp_contig_cnts = smalls[1]

ref_contig_dict = Counter(ref_contig_cnts)
ref_contigs = sorted(ref_contig_dict.keys(),reverse=True)
cp_ref = [sum(bigs[0])]

for contig in ref_contigs:
    cur = 0
    if len(cp_ref)>0:
        cur+=cp_ref[-1]
    cp_ref.append((cur+contig*ref_contig_dict[contig]))

del cp_ref[0]
    
cp_ref= list(map(lambda x:x/ref_tot,cp_ref))


exp_contig_dict = Counter(exp_contig_cnts)
exp_contigs = sorted(exp_contig_dict.keys(),reverse=True)
cp_exp = [sum(bigs[1])]
for contig in exp_contigs:
    cur = 0
    if len(cp_exp)>0:
        cur+=cp_exp[-1]
    cp_exp.append((cur+contig*exp_contig_dict[contig]))

del cp_exp[0]

cp_exp= list(map(lambda x:x/exp_tot,cp_exp))
    
#ref_contig_cnts.sort(reverse=True)
#exp_contig_cnts.sort(reverse=True)

#cp_ref = np.cumsum(ref_contig_cnts)/sum(ref_contig_cnts)
#cp_exp = np.cumsum(exp_contig_cnts)/sum(exp_contig_cnts)

x_ref = np.linspace(min(ref_contigs),max(ref_contigs),100)
x_exp = np.linspace(min(exp_contigs),max(exp_contigs),100)

xs_ref = np.array(ref_contigs[::-1])
ys_ref = np.array(cp_ref[::-1])
spl_ref = UnivariateSpline(xs_ref,ys_ref,k=1)

xs_exp = np.array(exp_contigs[::-1])
ys_exp = np.array(cp_exp[::-1])
spl_exp = UnivariateSpline(xs_exp,ys_exp,k=3)

f, ax = plt.subplots(1, 1, figsize=(10, 8))

ax.scatter(ref_contigs,cp_ref,c=colors[0],alpha = 0.7,s=3)
ax.scatter(exp_contigs,cp_exp,c=colors[1],alpha = 0.7,s=3)

ax.plot(x_ref,spl_ref(x_ref),color = colors[0],lw=3,alpha=0.3,linestyle = "--",label=leg_labels[0]+" Smoothed")
ax.plot(x_exp,spl_exp(x_exp),color = colors[1],lw=3,alpha=0.3,linestyle = "--",label=leg_labels[1]+" Smoothed")

ax.plot(ref_contigs,cp_ref,c=colors[0],alpha = 0.85,label=leg_labels[0])
ax.plot(exp_contigs,cp_exp,c=colors[1],alpha = 0.85,label=leg_labels[1])




ax.set_xscale('log')
ax.invert_xaxis()

ax.legend(prop={'size': 10})


plt.xlabel("Contig size",fontsize=12)
plt.ylabel("Proportion of genome",fontsize=12)
plt.title("cdf of proportion of genome filled vs smallest contig added\nOnly Contigs < 50KB",fontsize=15)
plt.show()
plt.close()



