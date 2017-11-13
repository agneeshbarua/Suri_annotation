#!/usr/bin/env python3

suri_file = open('surinamensis_nucleotide BM-83917.fasta.transdecoder.fasta')
file_with_adaptor = open('adaptors')

from Bio import SeqIO

adaptor_index_len = {}
adaptor_list = []

def trim_adaptors(seq_file, adaptor):
	len_adaptor = len(adaptor)
	for seq in seq_file:
		index = seq.seq.find(adaptor)
		if index == -1:
			#does not find adaptor, so wont trim
			yield seq
		else:
			#trim off the adaptor
			yield seq[index+len_adaptor:]
			#adaptor_index_len[index] = len_adaptor




for ad_seq in file_with_adaptor:
	ad_seq = ad_seq.rstrip()
	if not ad_seq.startswith('>'):
		adaptor_list.append(ad_seq)

print(adaptor_list)



original_reads = SeqIO.parse(suri_file, 'fasta')
#Need to find a way to loop this for each_seq in adaptor_list
trimmed_reads = trim_adaptors(original_reads, adaptor_list[0])
trimmed_reads2 = trim_adaptors(trimmed_reads, adaptor_list[1])
trimmed_reads3 = trim_adaptors(trimmed_reads2, adaptor_list[2])
trimmed_reads4 = trim_adaptors(trimmed_reads3, adaptor_list[3])
trimmed_reads5 = trim_adaptors(trimmed_reads4, adaptor_list[4])
trimmed_reads6 = trim_adaptors(trimmed_reads5, adaptor_list[5])
trimmed_reads7 = trim_adaptors(trimmed_reads6, adaptor_list[6])
trimmed_reads8 = trim_adaptors(trimmed_reads7, adaptor_list[7])
trimmed_reads9 = trim_adaptors(trimmed_reads8, adaptor_list[8])
trimmed_reads10 = trim_adaptors(trimmed_reads9, adaptor_list[9])
trimmed_reads11 = trim_adaptors(trimmed_reads10, adaptor_list[10])
trimmed_reads12 = trim_adaptors(trimmed_reads11, adaptor_list[11])
trimmed_reads13 = trim_adaptors(trimmed_reads12, adaptor_list[12])
trimmed_reads14 = trim_adaptors(trimmed_reads13, adaptor_list[13])
trimmed_reads15 = trim_adaptors(trimmed_reads14, adaptor_list[14])
trimmed_reads16 = trim_adaptors(trimmed_reads15, adaptor_list[15])
trimmed_reads17 = trim_adaptors(trimmed_reads16, adaptor_list[16])
trimmed_reads18 = trim_adaptors(trimmed_reads17, adaptor_list[17])
trimmed_reads19 = trim_adaptors(trimmed_reads18, adaptor_list[18])
trimmed_reads20 = trim_adaptors(trimmed_reads19, adaptor_list[19])
trimmed_reads21 = trim_adaptors(trimmed_reads20, adaptor_list[20])
trimmed_reads22 = trim_adaptors(trimmed_reads21, adaptor_list[21])
trimmed_reads23 = trim_adaptors(trimmed_reads22, adaptor_list[22])
trimmed_reads24 = trim_adaptors(trimmed_reads23, adaptor_list[23])
trimmed_reads25 = trim_adaptors(trimmed_reads24, adaptor_list[24])
trimmed_reads26 = trim_adaptors(trimmed_reads25, adaptor_list[25])

#Final output of trimmed RNAseq in fasta format
trimmed_file = SeqIO.write(trimmed_reads26, 'trimmed_suri_RNAseq.fa', 'fasta')






