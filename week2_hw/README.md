This homework assignment asks us to find a fasta file of an entire organisms genome sequence (less than 50MB)

I chose the genome sequence for psilocybe cyanescens (magic mushrooms)

We were to write a python program to read through the fasta file and output how many lines and headers the fasta file has.

We can verify the output of our program is correct by running the following commands to count number of lines:

wc -l <fasta_file_name>

And headers:

grep -c "^>" <fasta_file_name>

I ran my script & got a total of:

599274 lines

3976 headers

To use the python script:
 - I wrote it to work in python27 (most Linux distributions come with this version by default)
 - execute the script through python[27] count_lines_and_headers.py
   - can either send the name of the fasta file to the program through standard in, or leave the field blank & the script would default to look for and count the number of lines and headers located in this directory


