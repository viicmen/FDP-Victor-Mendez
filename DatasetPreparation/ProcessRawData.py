#
# FUNCTIONS FOR PROCESSING CISBP DATA
#

import os

def CisbpToMeme(input_dir, output_dir):
'''
Converts all motif files from the CIS-BP format in a directory to the MEME motif format and saves them in a new directory. 
Each output file includes a standard MEME header and the corresponding position weight matrix (PWM).
'''
    motifs = os.listdir(input_dir)

    for m in motifs:

        with open(input_dir+m, 'r') as f:
            matrix = []
            f.readline()
            
            for line in f.readlines():
                matrix.append(line.split()[1:])
        
        m = m[:-4]
        
        with open(output_dir + m + '.meme', 'w') as f:
        
            f.write('MEME version 4\n\n')
            f.write('ALPHABET= ACGT\n\n')
            f.write('Background letter frequencies\n')
            f.write('A 0.25 C 0.25 G 0.25 T 0.25\n\n')
            f.write(f'MOTIF {m}\n\n')
            f.write(f'letter-probability matrix: w= {str(len(matrix))}\n')
        
            for pos in matrix:
                f.write(' '.join(pos)+'\n')
            f.write('URL')

def WriteTfToM_Cisbp(out_file, input_file = 'TF_information.txt'):
'''
Parses a CIS-BP TF information file and creates a mapping from transcription factors (TFs) to their associated motif IDs. 
This mapping is saved to an output file and also returned as a dictionary.
'''
    dic = {}

    with open(input_file, 'r') as f:

        with open(out_file, 'w') as tf_to_m:

            f.readline()

            for line in f.readlines():

                line = line.split('\t')
                motifs = ''.join(line[6:]).split(',')
                motifs = motifs[:-1]
                ID = line [0]
                tf_to_m.write(ID + ' ' + ' '.join(motifs) + '\n')
                dic[ID] = motifs

    return dic

def Criteria(from_pos, to_pos):
'''
Checks whether a set of domains are located closely together by comparing the start and end positions of consecutive domains. 
Returns True if all are ≤40 amino acids apart, otherwise False.
'''
    for i in range(1,len(from_pos)):

        if int(from_pos[i]) - int(to_pos[i-1])> 40:
            return False

    return True

def CheckDomains(output_file, input_file = 'prot_seq.txt'):
'''
Reads protein sequence data with annotated domains and filters out TFs whose domains are too far apart. 
Writes the valid TFs, their domains, and sequences to an output file.
'''
    with open(input_file, 'r') as f, open(output_file, 'w') as out_file:

        for line in f.readlines():

            line = line.split()
            tf = line[0]
            domains = line[5].split(',')
            from_pos = line[6].split(',')
            to_pos = line[7].split(',')
            seq = line[9]

            if criteria(from_pos, to_pos):
                domains = ','.join(set(domains))
                out_file.write(tf + ' ' + domains + ' ' + seq + '\n')
    

def WriteFasta_Domains_Cisbp(tf_to_m, input_file, output_fasta, output_domains):
'''
Creates two output files from CIS-BP data: a FASTA file with sequences for each TF motif, and a file listing the corresponding domain annotations. 
It uses a TF-to-motif dictionary to resolve identifiers.
'''
    dic = read_tf_to_m(tf_to_m)

    with open(input_file, 'r') as f, open(output_fasta, 'w') as fasta, open(output_domains, 'w') as domains:

        for line in f.readlines():

            line = line.split()
            tf = line[0]
            domain = line[1]
            seq = line[2]

            if tf in dic:
                fasta.write(f'> {dic[tf]}\n{seq}\n')
                domains.write(f'{dic[tf]\t{domain}\n')

#
# FUNCTIONS FOR PROCESSING JASPAR DATA
#

import json

def WriteFasta_Jaspar(jaspar_json, out_fasta):
 '''
 This function creates a FASTA file from a JASPAR-formatted JSON file. 
 For each TF in the JSON, it writes one FASTA entry per PWM ID linked to that TF, using the full protein sequence.
 '''
    with open(jaspar_json) as f:
         dic = json.load(f)
         
    out_fasta = open(out_fasta, 'w')
    
    for tf in dic:
	       for m in dic[tf][0]:
              out_fasta.write(f'>{m}\n{dic[tf][1]}\n')

    out_fasta.close()
    
def WriteFastaDBD_TfToM_Jaspar(jaspar_json, dbds_path, pwm_dir, out_fasta, out_file):
'''
This function extracts domain sequences from full TF sequences using domain start and end positions, and saves them in a FASTA file. 
It also creates a file mapping TFs to their associated PWMs.
'''
    with open(jaspar_json) as f:
         dic = json.load(f)
         
    dbds = open(dbds_path, 'r')
    out_fasta = open(out_fasta, 'w')
    
    tftom = {}
    dbds.readline()
    pwm_dir = os.listdir(pwm_dir)

    for line in dbds.readlines():

        line = line.split()
        TF = line[0]
        PWM = line[1]
        FROM = int(line[2])
        TO = int(line[3])

        if PWM: 
        	out_fasta.write(f'>{PWM}\n{dic[TF][1][FROM:TO+1]}\n')
        	tftom[TF] = tftom.get(TF, []) + [PWM]

    dbds.close()
    out_fasta.close()
    
    with open(out_file, 'w') as f:
    		for tf in tftom:
    		    f.write(f'{tf}\t{'\t'.join(tftom[tf])}\n')

#
# FUNCTIONS FOR PROCESSING HOCOMOCO DATA
#

import json

def WriteFasta_Pwms_Hocomoco(json_hocomoco, uniprot_fasta, out_dir_pwms, out_fasta):
'''
This function reads motif data from a HOCOMOCO JSONL file and writes each PWM in MEME format to a specified directory. 
It also extracts protein sequences for each transcription factor from a UniProt FASTA file and writes them to a new FASTA file, linking each TF to its motif.
'''
	header = 'MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n'

	tfs = {}

	with open(json_hocomoco, "r") as f:

		for line in f:

		    data = json.loads(line)
		    motif = data['name']
		    l = data['length']
		    species = 'HUMAN' if 'HUMAN' in data['masterlist_info']['species'] else 'MOUSE'
		    tf = data['masterlist_info']['species'][species]['uniprot_ac']

		    tfs[tf] = tfs.get(tf, []) + [motif]

		    pwm = open(f'{out_dir_pwms}/{motif}.meme', 'w')
		    pwm.write(header)
		    pwm.write(f'MOTIF motif\n\n')
		    pwm.write(f'letter-probability matrix: alength= 4 w= {str(data['length'])}\n')
		    for pos in data['pfm']:

		        pos = [str(p) for p in pos]
		        pwm.write(f'{' '.join(pos)}\n')
		    pwm.write('URL')
	
    uni = open(uniprot_fasta, 'r')
    fasta = open(out_fasta, 'w')
    ok = False

    tf = fasta.readline().split('|')[1]
    if tf in tfs:
        new.write('> '+tfs[tf]+'\n')
        ok = True
    seq = ''

    for line in fasta.readlines():

        if line[0] == '>':
            if tf in tfs:
                new.write(seq[tfs[tf]['f']:tfs[tf]['t']+1]+'\n')
            tf = line.split('|')[1]
            if tf in tfs:

                new.write('>'+tfs[tf]['m']+'\n')
                ok = True
                seq = ''
            else:

                ok = False

        elif ok:

            seq+=line.strip()

    if tf in tfs:
        new.write(seq[tfs[tf]['f']:tfs[tf]['t']+1]+'\n')
    new.close()
    dbds.close()
    fasta.close()

def WriteFastaDBD_TfToM_Hocomoco(dbds_file, fasta_file, out_fasta, out_file):
'''
This function reads domain coordinates from a file and extracts the DNA-binding domain (DBD) region from protein sequences in a FASTA file with full TF sequences. 
It writes the DBD sequences to a FASTA file and generates a mapping file linking each transcription factor to its motif.
'''
    dbds = open(dbds_file,'r')
    f = open(out_file, 'w')
    tfs={}
    for line in dbds.readlines():
        line = line.split()
        TF = line[0]
        PWM = line[1]
        FROM = int(line[2])
        TO = int(line[3])
        tfs[TF]={'m':PWM, 'f':FROM,'t':TO}
        f.write(f'{TF}\t{'\t'.join(PWM)}\n')
    
    f.close()
    fasta = open(fasta_file, 'r')
    new = open(out_fasta, 'w')
    ok = False

    tf = fasta.readline().split('|')[1]
    if tf in tfs:
        new.write('> '+tfs[tf]+'\n')
        ok = True
    seq = ''

    for line in fasta.readlines():

        if line[0] == '>':
            if tf in tfs:
                new.write(seq[tfs[tf]['f']:tfs[tf]['t']+1]+'\n')
            tf = line.split('|')[1]
            if tf in tfs:

                new.write('>'+tfs[tf]['m']+'\n')
                ok = True
                seq = ''
            else:

                ok = False

        elif ok:

            seq+=line.strip()

    if tf in tfs:
        new.write(seq[tfs[tf]['f']:tfs[tf]['t']+1]+'\n')
    new.close()
    dbds.close()
    fasta.close()

#
# FUNCTIONS FOR CLASSIFYING TF FAMILIES WITH HHMER & PFAM
#

import pandas as pd
import os

def FindDomains(fasta_file, pfam_dbds, out_file):
'''
This function scans DNA-binding domain sequences using HMMER and a Pfam profile HMM database when executing "hmmer.sh". 
It parses the output and generates a TSV file listing the predicted domains for each sequence in the input FASTA files.
'''
    os.system(f'hmmscan --domtblout output -o dummy -E 0.001 {pfam_dbds} {fasta_file}')

    data = []

    with open('output', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split(maxsplit=22)
            if len(parts) == 23:
                data.append(parts)

    df = pd.DataFrame(data, columns = columns)

    grouped = df.groupby('query_name')['target_name'].unique().reset_index()
    results = pd.DataFrame({
    'seq_id': grouped['query_name'],
    'domains': grouped['target_name'].apply(lambda x: ', '.join(x))
})

    results.to_csv(out_file, sep='\t', index=False, na_rep='')


def ClassifyDomains(dbs):
'''
This function reads the domain assignments from the TSV files generated by FindDomains and organizes motifs into folders based on their domain type. 
It creates separate FASTA and MEME-formatted PWM files per domain, grouping and duplicating data accordingly for downstream analysis.
'''
    for db in dbs:

        print(db.split('/')[-2])
        dom_dic={}
        fasta = db+'filtered_sequences.fasta'

        with open(db+'sequence_domains.tsv', 'r') as f:
            header = f.readline()
            for line in f.readlines():

                motif, domain = line.split()
                dom_dic[domain] = dom_dic.get(domain, []) + [motif]

        for domain in dom_dic:

            print(domain,'\t', len(dom_dic[domain]))
            seq_path = db+'dbds/'+domain
            pwm_path = db+'pwms/'+domain

            if os.path.exists(pwm_path):

                os.system('rm -r '+seq_path)
                os.system('rm -r '+pwm_path)

            os.makedirs(seq_path)
            os.makedirs(pwm_path)

            for motif in dom_dic[domain]:

                os.system('cp '+db+'pwms/'+motif+'.meme '+pwm_path)
                os.system('cat '+db+'dbds/'+motif+'.fa >> '+seq_path+'/'+domain+'.fa')
                os.system('cat '+db+'dbds/'+motif+'.fa > '+seq_path+'/'+motif+'.fa')

        total = sum(len(dom_dic[domain]) for domain in dom_dic)
        print('total db\t',total)

