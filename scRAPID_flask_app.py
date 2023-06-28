from flask import Flask, request, Response
import csv
import sqlite3
from io import StringIO

app = Flask(__name__)

@app.route("/")
def hello():
    print("scrapid_db")
    return "scRAPID test"

# Please note that you'll need to replace the /download endpoint with the appropriate route in your Flask application.
@app.route('/database', methods=['GET'])
def download_csv():
    # Get query parameters from the request
    protein_query = request.args.get('protein')
    rna_query = request.args.get('rna')

    # Perform the query on the database to retrieve the desired data subset
    data_subset = query_database(protein_query, rna_query)

    # Generate the CSV file in-memory
    csv_data = generate_csv(data_subset)

    # Create a response with the CSV data
    response = Response(csv_data, mimetype='text/csv')
    response.headers.set('Content-Disposition', 'attachment', filename='catRAPID_table_subset.csv')

    return response

def query_database(protein_query, rna_query):
    # Implement the database query based on the provided protein and RNA queries
    # Modify the query based on the new query options for proteins and RNAs
    query = """
    SELECT * FROM catrapid
    WHERE (Uniprot_ID = ? OR Entry_Name = ?)
    AND (Ensembl_Gene_ID = ? OR Ensembl_Transcript_ID = ? OR gene_name = ?)
    """

    # Create a connection to the SQLite database
    conn = sqlite3.connect('catrapid.db')
    cursor = conn.cursor()

    try:
        # Execute the query with the provided parameters
        cursor.execute(query, (protein_query, protein_query, rna_query, rna_query, rna_query))

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
    csv_data = StringIO()
    writer = csv.writer(csv_data)

    # Write the CSV header
    writer.writerow(['Uniprot_ID', 'Entry_Name', 'Protein_Frag', 'Ensembl_Gene_ID', 'Ensembl_Gene_ID_version', 'Ensembl_Transcript_ID', 'Ensembl_Transcript_ID_version', 'gene_name', 'RNA_Frag', 'score'])

    # Write the data rows
    for row in data_subset:
        writer.writerow(row)

    return csv_data.getvalue()
