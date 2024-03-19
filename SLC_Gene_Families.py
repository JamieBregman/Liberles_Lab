"""
Script for single linkage clustering (SLC)
"""

### IMPORTING MODULES ###
import sys
import csv
import pandas as pd

### PARSING INPUT CSV ###
#Creating a dictionary value for each row based on the threshold value
csv_file = sys.argv[1] #CSV input
hits_dict = {} #Empty dictionary to store hit pairs

#Opening CSV file and adding each row to a dictionary
#Format: {[index]:[gene1, gene2, %ID, %Ungapped]}
with open(csv_file, "r") as blast_hits:
    file = csv.reader(blast_hits, delimiter = ",")

    #Initializing counter (counter becomes dictionary key)
    counter = 0

    #Iterating over each row in the CSV to remove values below the threshold and creating the dictionary
    for row in file:
        if float(row[2]) >= float(sys.argv[2]) and float(row[3]) >= float(sys.argv[3]): #If both values are above the threshold
            hits_dict[counter] = [row[0], row[1], float(row[2]), float(row[3])] #Add values
            counter += 1 #Increase counter


### CREATING LISTS/DICTS NEEDED FOR SLC ###
#Creating a list of all unique genes
temp_gene_list = [] #Temporary gene list
#Adding the first 2 values for each dictionary key (the genes)
for key in hits_dict:
    temp_gene_list.extend(hits_dict[key][:2])

gene_list = list(set(temp_gene_list)) #Using set to create a list of alll unique genes in the MAFFT input file

index_dict = {index:item for index, item in enumerate(gene_list)} #creating master index dictionary (assigning each gene a unique index)
index_dict_reversed = {value:key for key, value in index_dict.items()} #creating reversed master index dictionary (used to create final gene list at the end of the file)

#Creating master family list
fam_list = list(index_dict.keys()) #List of all unique indexes for each gene


### MAIN LOOP ###
#For each gene pair:
for pair in hits_dict:
    val0 = index_dict_reversed[hits_dict[pair][0]] #val0 = unique index of gene0
    val1 = index_dict_reversed[hits_dict[pair][1]] #val1 = unique index of gene1

    #If gene0 and gene1 are the same
    if val0 == val1:
        continue
    else:
        #If gene0 index > gene1 index
        if val0 >= val1:
            temp_val0 = fam_list[val0] #temp_val0 = value that was in the index of gene0
            fam_list[val0] = fam_list[val1] #replace value in gene0 index with value in gene1 index
            fam_list = [fam_list[val1] if i == temp_val0 else i for i in fam_list] #replace all other occurrences of the value that was in index of gene0

        #If gene1 index > gene0 index
        else:
            temp_val1 = fam_list[val1]
            fam_list[val1] = fam_list[val0]
            fam_list = [fam_list[val0] if i == temp_val1 else i for i in fam_list]
            

### EDITING fam_list TO CREATE FINAL FAMILIES ###
final_dict = {val:[] for val in set(fam_list)} #dictionary of all unique items in the family list after looping

#for each item in the family list, append that item's index to it's corresponding key in the dictionary
for index,item in enumerate(fam_list):
    final_dict[item].append(index)

fam_list_of_lists = list(final_dict.values()) #creating a list of all of the sublists from the family list
fam_list_FINAL = [[index_dict.get(item, item) for item in sublist] for sublist in fam_list_of_lists] #replacing the indexes with the gene names
#print(fam_list_FINAL)

### WRITING FAMILIES TO CSV FILE WHERE EACH ROW IS A UNIQUE FAMILY ###
with open("SLC_OUTPUT.txt", "w") as output_file:
    for sublist in fam_list_FINAL:
        fam = " ".join(map(str, sublist)) + "\n"
        output_file.write(fam + "\n")