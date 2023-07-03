from flask import Flask, request, Response
import csv
import sqlite3
from io import StringIO

app = Flask(__name__)

@app.route('/database', methods=['GET'])
def download_csv():
    # Get query parameters from the request
    protein_query = request.args.get('protein')
    rna_query = request.args.get('rna')

    # Perform the query on the database to retrieve the desired data subset
    data_subset = query_database(protein_query, rna_query)

    # Generate the CSV file
    csv_file = generate_csv(data_subset)

    # Create a response with the CSV file
    response = Response(csv_file.getvalue(), mimetype='text/csv')
    response.headers.set('Content-Disposition', 'attachment', filename='catRAPID_table_subset.csv')

    return response

def query_database(protein_query, rna_query):
    # Create a connection to the SQLite database
    conn = sqlite3.connect('catrapid.db')
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
    csv_file = StringIO()
    writer = csv.writer(csv_file)

    # Write the CSV header
    writer.writerow(['Uniprot_ID', 'Entry_Name', 'Protein_Frag', 'Ensembl_Gene_ID', 'Ensembl_Gene_ID_version', 'Ensembl_Transcript_ID', 'Ensembl_Transcript_ID_version', 'gene_name', 'RNA_Frag', 'score'])

    # Write the data rows
    writer.writerows(data_subset)

    return csv_file
