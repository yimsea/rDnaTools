from pbcore.io.FastaIO import FastaReader, FastaWriter

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

def copy_fasta_list( sequence_list, output_file ):
    with FastaWriter( output_file ) as writer:
        with open( sequence_list ) as handle:
            for line in handle:
                sequence_file = line.strip()
                copy_fasta_sequences( sequence_file, writer )

if __name__ == '__main__':
    import sys

    sequence_list = sys.argv[1]
    output_file = sys.argv[2]

    copy_fasta_list( sequence_list, output_file )