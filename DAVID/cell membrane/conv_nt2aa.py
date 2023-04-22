from Bio.Seq import translate

fasta = open('cell_membrane_transcripts.fasta', 'r')

for line in fasta:
	if line.startswith('>'):
		print(line)
	else:
		print(translate(line))
