import sys
filename = "mytestfile.fastq"

fastq_file = open(filename)

lines = fastq_file.read().strip().split("\n")
print(lines[3])
scores = [ord(c)-33 for c in lines[3]]
print(scores)
probs = [10.0**(-1*q/10.0) for q in scores]
avg_prob_of_wrong = sum(probs)/len(probs)
print("Sequence of length: {}".format(len(probs)))
print("prob incorrect: {}".format(avg_prob_of_wrong))
print("prob correct: {}".format(1-(sum(probs)/len(probs))))

"""
for i,line in enumerate(sys.stdin):
    if i%4 == 3:
        avg_prob_of_wrong = sum([10.0**(-1*(ord(c)-33)/10.0) for c in line.strip()])/len(line.strip())
        print("\navg prob that seq {} is:\n   incorrrect: {}\n      correct: {}".format(i//4+1,avg_prob_of_wrong,1-avg_prob_of_wrong))
"""
