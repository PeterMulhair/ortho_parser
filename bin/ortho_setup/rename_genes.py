import os
import glob
from Bio import SeqIO
from subprocess import call as unix

'''
Script to rename primary transcript genes
to use as input for OrthoFinder
'''

os.makedirs('genes_renamed', exist_ok=True)

for fa in glob.glob('*fa'):
    spName = fa.split('.')[0]
    spName_shr = spName.split('_')[0][:3] + spName.split('_')[1][:4].title()
    print(spName_shr)
    with open(fa) as f, open('genes_renamed/' + spName_shr + '.fa', 'w') as outF:
        for record in SeqIO.parse(f, 'fasta'):
            header = record.description
            seq = str(record.seq)
            outF.write('>' + spName_shr + '_' + header + '\n' + seq + '\n')
