"""
Extract hit pairs from BLAST
Perform MAFFT pairwise alignment
Record %id & %ungapped for each hit pair
Save quality control values to file
"""

###IMPORTING MODULES###
import MySQLdb
import sys
import csv
import subprocess
import time
from Bio.Align.Applications import MafftCommandline

### ACCESSING MYSQL DB AND EXTRACTING SEQUENCES ###
#Connecting to species_protein_seqs database
mydb = MySQLdb.connect(
                host ='localhost',
                user='jbregman',
                passwd='jbregman',
                db='species_protein_seqs')

mycursor = mydb.cursor() #Initializing cursor
mycursor.execute("SELECT * FROM seqs;") #Acessing seqs table
mydb.commit()

table_rows = mycursor.fetchall() #Storing each row of seqs as a tuple
sql_dict = {item[1]:item[3] for item in table_rows} #Creating a dictionary from the SQL table

### ACCESSING BLAST OUTPUT CSV FILE AND EXTRACTING HIT PAIRS ###
csv_file = sys.argv[1] 
hits_dict = {} #Empty dictionary to store hit pairs

#Opening CSV file and adding the hit pairs to a dictionary
with open(csv_file, "r") as blast_hits:
    file = csv.reader(blast_hits, delimiter = ",")

    counter = 0
    for row in file:
        hits_dict[counter] = [row[0], row[1]]
        counter += 1


### MAIN LOOP###
output_dictionary = {} #Creating dicitionary used to add values to the CSV

### CREATING A FASTA FROM EACH HIT PAIR###
for pair in hits_dict:

    #Extracting hit pair IDs from the hits dictionary
    hit_id1 = hits_dict[pair][0]
    hit_id2 = hits_dict[pair][1]

    if hit_id1 == hit_id2:

        #Calculating the values
        percent_id = float(1.0)
        percent_ungapped = float(1.0)

        #Adding values to dictionary
        output_dictionary[pair] = [percent_id, percent_ungapped]

    else:
        #Extracting hit pair sequences from SQL dictionary
        hit_seq1 = sql_dict[hit_id1]
        hit_seq2 = sql_dict[hit_id2]

        #Creating fasta file
        temp_fasta = open("%s" % pair + ".fa", "w")
        temp_fasta.write(">" + str(hit_id1) + "\n")
        temp_fasta.write(str(hit_seq1) + "\n")
        temp_fasta.write(">" + str(hit_id2) + "\n")
        temp_fasta.write(str(hit_seq2))
                                                                                
        temp_fasta.close()
    
        ### RUNNING MAFFT AND PARSING OUTPUT ###
        mafft_path = "/usr/bin/mafft" #Path to mafft executable
        mafft_command = MafftCommandline(mafft_path, input = "%s" % pair + ".fa", amino = True) #Running mafft
        stdout, stderr = mafft_command() #Storing mafft output

        #Parsing mafft output
        mafft_output_split = stdout.split("\n")

        temp_ids = []
        temp_seqs = []
        current_id = None
        current_seq = ""
    
        for item in mafft_output_split:
            if item.startswith(">"):
                if current_id is not None:
                    temp_seqs.append(current_seq)
                    current_seq = ""
                current_id = item
                temp_ids.append(item)
            else:
                current_seq += item

        if current_seq:
            temp_seqs.append(current_seq)

        ### CALCULATING % VALUES FOR THE ALIGNED SEQUENCES ###  
        #Initalizing values
        match_count = 0
        mismatch_count = 0
        gapped_count = 0
    
        #Looping through the mafft aligned sequences and calculating %id and %ungapped
        for i in range(len(temp_seqs[0])):
            if temp_seqs[0][i] == "-":
                if temp_seqs[1][i] == "-":
                    continue
                else:
                    gapped_count += 1
            else:
                if temp_seqs[0][i] == temp_seqs[1][i]:
                    match_count += 1
                else:
                    if temp_seqs[1][i] == "-":
                        gapped_count += 1
                    else:
                        mismatch_count += 1

        #Calculating the values
        percent_id = float(match_count)/float(match_count + mismatch_count)
        percent_ungapped = float(match_count + mismatch_count)/float(match_count + mismatch_count + gapped_count)

        #Adding values to dictiionary
        output_dictionary[pair] = [percent_id, percent_ungapped]

        #Removing the created fasta file
        subprocess.call("rm *.fa", shell = True)


### ADDING THE %ID AND %UNGAPPED VALUES TO A CSV FILE ###
#Reading in the CSV data
with open(sys.argv[1], "r") as csv_input:
    csv_read = csv.reader(csv_input)
    csv_data = list(csv_read)

#Formatting the %id and %ungapped values
for i, row in enumerate(csv_data):
    if i in output_dictionary:
        value1, value2 = output_dictionary[i]
        row.extend([value1, value2])

#Writing the values to the CSV
with open("NEW_%s" % sys.argv[1], "w") as csv_output:
    writer = csv.writer(csv_output)
    writer.writerows(csv_data)
