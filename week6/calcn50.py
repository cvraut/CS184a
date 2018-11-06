"""
generate lengths file:
    bioawk -cfastx '{print $name,length($seq)}' fastq_name > lengths.txt

    this py script expects input to be from the format of the above basnh cmd
"""

import sys

contigs = []
tot_sum = 0
for line in sys.stdin:
    contig_name,length = line.split()
    length = int(length)
    contigs.append((length,contig_name))
    tot_sum += length

contigs.sort(reverse = True)
cur_tot = 0

ng50 = 150e6/2
for l,c in contigs:
    cur_tot+=l
    ng50-=l
    if cur_tot>tot_sum/2:
        print("N50 is: {}".format(l))
        cur_tot = 0
    if ng50 < 0:
        print("NG50 is: {}".format(l))
        ng50 = 1e15

hacked_contigs_10000 = filter(lambda x: x[0]>10000, contigs)
hacked_contigs_1000 = filter(lambda x: x[0]>1000, contigs)
hacked_contigs_mil = filter(lambda x: x[0]>1000000, contigs)
#print(len(contigs),len(hacked_contigs_1000),len(hacked_contigs_10000))
#print(tot_sum,sum(map(lambda x:x[0],hacked_contigs_1000)),sum(map(lambda x:x[0],hacked_contigs_10000)))

print("\nFiltering contigs greater than 1000:")
tot_sum = sum(map(lambda x:x[0],hacked_contigs_10000))
cur_tot = 0

ng50 = 150e6/2
for l,c in hacked_contigs_1000:
    cur_tot+=l
    ng50-=l
    if cur_tot>tot_sum/2:
        print("N50 is: {}".format(l))
        cur_tot = 0
    if ng50 < 0:
        print("NG50 is: {}".format(l))
        ng50 = 1e15

print("\nFiltering contigs greater than 10000:")
tot_sum = sum(map(lambda x:x[0],hacked_contigs_1000))

cur_tot = 0

ng50 = 150e6/2
for l,c in hacked_contigs_10000:
    cur_tot+=l
    ng50-=l
    if cur_tot>tot_sum/2:
        print("N50 is: {}".format(l))
        cur_tot = 0
    if ng50 < 0:
        print("NG50 is: {}".format(l))
        ng50 = 1e15

print("\nFiltering contigs greater than 1e6:")
tot_sum = sum(map(lambda x:x[0],hacked_contigs_mil))

cur_tot = 0

ng50 = 150e6/2
for l,c in hacked_contigs_mil:
    cur_tot+=l
    ng50-=l
    if cur_tot>tot_sum/2:
        print("N50 is: {}".format(l))
        cur_tot = 0
    if ng50 < 0:
        print("NG50 is: {}".format(l))
        ng50 = 1e15


