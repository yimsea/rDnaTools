from pbcore.io.FastaIO import FastaReader

def fasta_count( fasta_file ):
    count = 0
    try:
        for record in FastaReader( fasta_file ):
            if len(record.sequence) > 0:
                count += 1
    except:
        pass
    return count

def copy_fasta_sequences( fasta_file, fasta_writer ):
    for fasta_record in FastaReader( fasta_file ):
        fasta_writer.writeRecord( fasta_record )
