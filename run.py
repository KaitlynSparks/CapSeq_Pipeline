if __name__ ==  '__main__':
    import uuid
    import os
    import glob
    import pybedtools
    import useful_functions
    import pandas as pd
    import subprocess
    from Bio import SeqIO
    from datetime import datetime

    #job parameters
    filter_for_capseq_results = True
    filter_for_genes = True
    filter_for_transcripts = False
    clear_temp = False

    data_directory = "data/"
    input_files_directory = "input/"
    output_directory = "jobs/"
    
    
    #create random job id
    job_id = str(uuid.uuid4())
    print("Started job: " + job_id)
    
    #create directories
    temp_folder = output_directory + job_id + "/temp_files/"
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    results_folder = output_directory + job_id + "/results/"
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    #filtering capseq results
    if filter_for_capseq_results:
        with open('input/CaptureSEQ_results_of_interest.txt') as f:
            lines = f.read().splitlines()
        results_bedfile_name_list = lines
    else: #needs testing
        results_bedfile_name_list = glob.glob(data_directory + "CaptureSEQ_results/"+ '*.bed')
        results_bedfile_name_list = map(list.split('/')[-1], results_bedfile_name_list)

    #iterate over each capseq results file
    for results_bedfile_name in results_bedfile_name_list:
        results_bedtool = pybedtools.BedTool(data_directory + "CaptureSEQ_results/" + results_bedfile_name)
        gene_masked_bedfile_name = temp_folder + 'gene_masked_' + results_bedfile_name
        results_bedtool = results_bedtool.sort(genome="hg38")
        
        #filter for genes, if applicable
        if filter_for_genes:
            with open(input_files_directory + 'genes_of_interest.txt', encoding="utf-8") as f:
                lines = f.read().splitlines()
            useful_functions.gene_filter(lines)
            gene_masking_bedtool = pybedtools.BedTool(temp_folder + 'gene_masking_bedfile.bed')
            results_bedtool = results_bedtool.intersect(gene_masking_bedtool, f=0.9, u=True, s=True)
        results_bedtool.saveas(gene_masked_bedfile_name)

        gene_masked_bed = pd.read_csv(gene_masked_bedfile_name, sep = '\t')

        #filter for transcripts, if applicable
        transcript_and_gene_masked_bed = gene_masked_bed
        if filter_for_transcripts:
            with open(input_files_directory + 'transcripts_of_interest.txt', encoding="utf-8") as f:
                lines = f.read().splitlines()
            transcript_and_gene_masked_bed = useful_functions.transcript_filter(transcript_and_gene_masked_bed, lines)
        transcript_and_gene_masked_bed_file_name = temp_folder + 'transcript_and_' + gene_masked_bedfile_name.split('/')[-1]
        transcript_and_gene_masked_bed.to_csv(transcript_and_gene_masked_bed_file_name, sep='\t', index=False, header=False)

        #extract DNA sequences
        final_bedtool = pybedtools.BedTool(transcript_and_gene_masked_bed_file_name)
        dna_fasta_file_name = temp_folder + 'dna_sequence_' + results_bedfile_name[0:-4] + '.fa'
        dna_sequences = final_bedtool.sequence(fi= 'data/hg38.fa', s=True, split=True, fo=dna_fasta_file_name)

        #get amino acid sequences (move to functions file at some point)
        dna_sequences = list(SeqIO.parse(dna_fasta_file_name, "fasta"))
        for i in range(len(dna_sequences)):
            all_reading_frames_aa_sequence_fasta_file = temp_folder + 'all_reading_frames_' + str(i) + '_' + results_bedfile_name[0:-4] + '.fa'
            cmd = 'curl -s -d "dna_sequence=' + str(dna_sequences[i].seq) + '&output_format=fasta" https://web.expasy.org/cgi-bin/translate/dna2aa.cgi > ' + all_reading_frames_aa_sequence_fasta_file
            subprocess.run(cmd, shell=True)
            aa_sequences = list(SeqIO.parse(all_reading_frames_aa_sequence_fasta_file, "fasta"))
            if transcript_and_gene_masked_bed.iloc[i,5] == "+": 
                del aa_sequences[0:3]
            elif transcript_and_gene_masked_bed.iloc[i,5] == "-": #neg strand not functional (to do with index?)
                del aa_sequences[4:7]
            else:
                print("No strand information found for " + transcript_and_gene_masked_bed.iloc[i,3] + " in " + all_reading_frames_aa_sequence_fasta_file )
            longest_length = 0
            longest_sequence = None
            for j in aa_sequences:
                split_aa_sequences = str(j.seq).split("-")
                chopped_aa_sequences = map(useful_functions.chop_sequence, split_aa_sequences)
                sorted_aa_sequences= sorted(chopped_aa_sequences, key=len)
                if len(sorted_aa_sequences[-1]) > longest_length:
                    longest_sequence = sorted_aa_sequences[-1]
                    longest_length = len(sorted_aa_sequences[-1])
            assert(longest_sequence[0] == "M")
            with open(results_folder + str(transcript_and_gene_masked_bed.iloc[i,3].split(";")[1]) + "_" + results_bedfile_name[0:-4] + ".txt", "w") as text_file:
                text_file.write(longest_sequence)

    #output amino acid sequences to combined fasta file (move to functions file at some point)
    results_bedfile_name_list = glob.glob(results_folder + '*.txt')
    sequence_list = []
    name_list = []
    for i in results_bedfile_name_list:
        with open(i,"r") as input_file:
            sequence = input_file.read()
            name = str(i.split("/",)[2])[0:-4]
            if sequence in sequence_list:
                name_list[sequence_list.index(sequence)] = name_list[sequence_list.index(sequence)] + "/" +name
            else:
                sequence_list.append(sequence)
                name_list.append(name)
    
    assert(len(sequence_list) == len(name_list))

    with open(results_folder + "combined_aa_sequences.fa","w") as output_file:
        for k in range(len(sequence_list)):
            if len(sequence_list[k])>400:
                output_file.write(">CU|" + name_list[k] + "\n")
                output_file.write(sequence_list[k] + "\n")
        
    #write parameters used in job to txt file
    with open("parameter_file.txt", "w") as parameter_file:
        parameter_file.write("submitted: " + str(datetime.now()))
        parameter_file.write("filter_for_capseq_results = " + str(filter_for_capseq_results))
        parameter_file.write("filter_for_genes = " + str(filter_for_genes))
        parameter_file.write("filter_for_transcripts = " + str(filter_for_transcripts))
        parameter_file.write("data_directory = " + data_directory)
        parameter_file.write("input_files_directory = " + input_files_directory)
        parameter_file.write("output_directory = " + output_directory)

    if clear_temp:
        files = glob.glob(temp_folder + '*')
        for f in files:
            os.remove(f)

