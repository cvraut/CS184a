fasta_file_name = "GCA_002938375.1_Psicy2_genomic.fna"

user_in = raw_input("Please enter the name of the fasta file you want to count (empty for GCA_002938375.1_Psicy2_genomic.fna): ")

if not user_in == "":
    fasta_file_name = user_in

fasta_file = open(fasta_file_name)

lines = fasta_file.readlines()
headers = sum([int(l.strip()[0]=='>') for l in lines])

print("The number of lines is: {}\nThe number of headers is: {}".format(len(lines),headers))

fasta_file.close()
