from flask import Flask, request, Response
import os
import csv
import sqlite3
import pandas as pd
from io import StringIO

app = Flask(__name__)

data_folder = "/Users/alexandrosarmaos/projects_local/scrapid/"
data_folder = "/var/www/html/scRAPID/database"


@app.route('/', methods=['GET'])
def hello():
    return "Hello"

@app.route('/database', methods=['GET'])
def download_csv():
    # Get query parameters from the request
    if 'protein' in request.args:
        protein_query = request.args.get('protein')
    else:
        protein_query = None
    
    if 'rna' in request.args:
        rna_query = request.args.get('rna')
    else:
        rna_query = None

    protein_uniprot_ids = None 

    if protein_query:
        protein_uniprot_ids = lookup_protein_query(protein_query)
    
    if protein_uniprot_ids:
        # Perform the query on the database to retrieve the desired data subset
        data_subset = query_database(protein_uniprot_ids, rna_query)
    else: 
        data_subset = query_database(protein_query, rna_query)
    
    # Generate the CSV file
    csv_file = generate_csv(data_subset)

    
    # Create a response with the CSV file
    response = Response(csv_file.getvalue(), mimetype='text/csv')
    response.headers.set('Content-Disposition', 'attachment', filename='catRAPID_table_subset.csv')

    return response

def lookup_protein_query(protein_query):
    lookup=pd.read_csv(os.path.join(data_folder, 'lookup_table_proteins.csv'))
    if protein_query in list(lookup.protein_name):
        return lookup[lookup.protein_name == protein_query]['accession_number'].values[0]
    else:
        return None
    
def query_database(protein_query, rna_query):
    # Create a connection to the SQLite database
    conn = sqlite3.connect(os.path.join(data_folder, 'catrapid.db'))
    cursor = conn.cursor()

    try:
        if protein_query and rna_query:
            query = """
            SELECT * FROM catrapid
            WHERE (Uniprot_ID = ? OR Entry_Name = ?)
            AND (Ensembl_Gene_ID = ? OR Ensembl_Transcript_ID = ? OR gene_name = ?)
            """
            cursor.execute(query, (protein_query, protein_query, rna_query, rna_query, rna_query))
        elif protein_query:
            # Single protein
            query = '''
            SELECT * FROM catrapid
            WHERE Uniprot_ID = ? OR Entry_Name = ?
            '''
            cursor.execute(query, (protein_query, protein_query))

        elif rna_query:
            # Single RNA
            query = '''
            SELECT * FROM catrapid
            WHERE Ensembl_Gene_ID = ? OR Ensembl_Transcript_ID = ? OR gene_name = ?
            '''
            cursor.execute(query, (rna_query, rna_query, rna_query))

        # Fetch all rows from the result set
        data_subset = cursor.fetchall()

        return data_subset

    except Exception as e:
        # Handle any errors that occurred during the execution
        print(f"Error executing database query: {e}")

    finally:
        # Close the cursor and the database connection
        cursor.close()
        conn.close()

def generate_csv(data_subset):
    lookup = pd.read_csv(os.path.join(data_folder,'lookup_table_proteins.csv'))
    csv_file = StringIO()
    writer = csv.writer(csv_file)

    # Write the CSV header
    writer.writerow(['Uniprot_ID', 'Entry_Name', 'Protein_Frag', 'Ensembl_Gene_ID', 'Ensembl_Gene_ID_version', 'Ensembl_Transcript_ID', 'Ensembl_Transcript_ID_version', 'gene_name', 'RNA_Frag', 'score', 'protein_name'])

    # Write the data rows
    for row in data_subset:
        # Retrieve the protein name based on the Uniprot ID
        protein_name=""
        if row[0] in list(lookup.accession_number):
            protein_name = lookup[lookup.accession_number==row[0]]['protein_name'].values[0]

        # Append the protein name to the row
        row_with_protein_name = list(row) + [protein_name]
        
        # Write the row to the CSV file
        writer.writerow(row_with_protein_name)

    return csv_file

#if __name__ == '__main__':
#    app.run()
