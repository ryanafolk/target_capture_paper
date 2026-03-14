#!/usr/bin/env python3

# NCBI SEQUENCE DOWNLOAD FOR SYMBIOSIS GENES
# This script downloads protein sequences of key symbiosis genes 
# from selected genera using NCBI Entrez. It ensures:
# - maximum sequences per gene (MAX_PER_GENE)
# - no duplicate sequences (by accession)
# - rescue search in all bacteria if none found in target genera

from Bio import Entrez, SeqIO
from time import sleep
from collections import defaultdict

# REQUIRED BY NCBI
# NCBI requires you to provide a valid email for queries
Entrez.email = "pdg99@msstate.edu"

# OUTPUT
# All downloaded sequences will be saved in this FASTA file
output_fasta = "symbiosis_all_genes_max1000.faa"

# FULL LIST OF GENES TO DOWNLOAD
genes = [
    "fdxB","fdxN",
    "fixA","fixB","fixC","fixJ","fixK","fixL","fixN","fixO","fixP","fixQ","fixU","fixX",
    "Irr","nfeD",
    "nifA","nifB","nifD","nifE","nifH","nifJ","nifK","nifM","nifN","nifQ","nifS","nifT","nifU","nifV","nifW","nifX","nifY","nifZ",
    "nodA","nodB","nodC","nodD","nodE","nodF","nodH","nodI","nodJ","nodL","nodP","nodQ","nodS","nodT","nodU","nodX","nodZ",
    "noeA","noeB","noeE","noeI","noeJ",
    "nolB","nolF","nolG","nolL","nolU","nolX",
    "nopE","nopL","nopM","nopP","nopT",
    "rhcN","rhcQ","rhcR","rhcS","rhcT",
    "rpoN"
]

# FULL LIST OF TARGET GENERA
genera = [
    "Methylobacterium",
    "Burkholderia",
    "Sinorhizobium",
    "Rhizobium",
    "Cupriavidus",
    "Azorhizobium",
    "Frankia",
    "Mesorhizobium",
    "Bradyrhizobium",
    "Ensifer",
    "Bosea",
    "Ochrobactrum",
    "Devosia",
    "Microvirga",
    "Aminobacter",
    "Phyllobacterium",
    "Neorhizobium",
    "Pararhizobium",
    "Shinella"
]

# PARAMETERS
MAX_PER_GENE = 1000   # maximum sequences per gene
BATCH_SIZE = 100      # sequences fetched per query
SLEEP_TIME = 1        # pause between queries (seconds)
MAX_RETRIES = 5       # retry attempts in case of failure

# DATA STRUCTURES TO TRACK SEQUENCES
seen_accessions = set()        # avoid duplicates
records_out = []               # store SeqIO records for output
gene_counts = defaultdict(int) # count sequences per gene

# FUNCTION TO DOWNLOAD SEQUENCES FOR ONE GENE AND GENUS
def download_sequences(gene, genus, retmax):
    """
    Download sequences from NCBI Protein database for a given gene and genus.
    Avoids duplicate accessions and updates records_out and seen_accessions.
    """
    for attempt in range(MAX_RETRIES):
        try:
            # Search NCBI Protein database for gene in genus
            handle = Entrez.esearch(
                db="protein",
                term=f'{gene}[Gene Name] AND "{genus}"[Organism]',
                retmax=retmax
            )
            record = Entrez.read(handle)
            handle.close()

            # If no sequences found, return 0
            if not record["IdList"]:
                return 0

            # Fetch sequences in FASTA format
            fetch = Entrez.efetch(
                db="protein",
                id=",".join(record["IdList"]),
                rettype="fasta",
                retmode="text"
            )

            count = 0
            # Parse FASTA and add to records_out
            for seq in SeqIO.parse(fetch, "fasta"):
                acc = seq.id.split("|")[0]
                if acc in seen_accessions:
                    continue
                seen_accessions.add(acc)
                count += 1
                # Clean header: gene|genus|accession
                seq.id = f"{gene}|{genus}|{acc}"
                seq.description = ""
                records_out.append(seq)

            fetch.close()
            sleep(SLEEP_TIME)
            return count

        except Exception as e:
            print(f"Retry {attempt+1}/{MAX_RETRIES} failed for {gene} in {genus}: {e}")
            sleep(SLEEP_TIME * 5)

    print(f"Failed after {MAX_RETRIES} retries: {gene} in {genus}")
    return 0

# MAIN LOOP: DOWNLOAD SEQUENCES FOR ALL GENES
for gene in genes:
    print(f"Downloading gene: {gene}")

    # Download sequences from target genera
    for genus in genera:
        remaining = MAX_PER_GENE - gene_counts[gene]
        if remaining <= 0:
            break
        count = download_sequences(gene, genus, min(BATCH_SIZE, remaining))
        if count > 0:
            gene_counts[gene] += count
            print(f"{count} sequences from {genus} added (total for {gene}: {gene_counts[gene]})")

    # Rescue mode: if no sequences found in target genera, search all bacteria
    if gene_counts[gene] == 0:
        print(f"Rescue search for {gene} in all bacteria...")
        count = download_sequences(gene, "bacteria", MAX_PER_GENE)
        if count > 0:
            gene_counts[gene] += count
            print(f"{count} sequences rescued for {gene} (total: {gene_counts[gene]})")
        else:
            print(f"No sequences found for {gene} in all bacteria.")

# WRITE OUTPUT FASTA
SeqIO.write(records_out, output_fasta, "fasta")

# FINAL SUMMARY
print("\nDOWNLOAD COMPLETE")
print(f"Total sequences downloaded: {len(records_out)}")
print(f"Output file: {output_fasta}\n")

print("Sequences per gene:")
for gene in genes:
    print(f"{gene}: {gene_counts[gene]}")