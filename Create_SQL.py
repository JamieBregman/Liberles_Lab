#Importing modules
from Bio import SeqIO
import MySQLdb
import re

#Connecting to mysql to access the database
mydb = MySQLdb.connect(
        host ='localhost',
        user='jbregman',
        passwd='jbregman',
        db='species_protein_seqs')

mycursor = mydb.cursor()
mydb.commit()

#Creating a table for the NCBI information
mycursor.execute("""CREATE TABLE seqs_new(
                unique_ID INT AUTO_INCREMENT PRIMARY KEY,
                sequence_ID VARCHAR(300),
                species_name VARCHAR(100),
                sequence TEXT);""")
                
mydb.commit()

#Adding fasta to this newly created table
with open("ALL_SPECIES.faa","r") as fasta:

    row_count = 0

    for record in SeqIO.parse(fasta, 'fasta'):
        sequence_id = record.id
        sequence = str(record.seq)

        match = re.findall(r'\[(.*?)\]', record.description)
        if match:
            species_name = match[-1]
        else:
            species_name = None
        
        mycursor.execute("""
        INSERT INTO seqs_new (sequence_id, species_name, sequence)
        VALUES (%s, %s, %s);""", (sequence_id, species_name, sequence))

        row_count += 1
        print(row_count)

mydb.commit()
mydb.close()
