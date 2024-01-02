
import os
import sys
import subprocess
from Bio import SeqIO
from collections import defaultdict

###################################################################################
# Autor: Jose Cleydson F Silva
# Script to compare large and complex genomes
# Versao: 1.2
###################################################################################

def Renome_seq(key,aux):
    aux.pop(key, None)
    return aux
    
##############################################################################
fasta = {}
with open(sys.argv[1], "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        fasta[record.id] = record.seq
##############################################################################
fastah = {}
with open(sys.argv[2], "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        fastah[record.id] = record.seq
###############################################################################

for key in sorted(fasta):
    print ('*************************** --- *************************** --- *************************** --- ***************************')
    print ('Query: '+ key+'\n') 
    seq = fasta[key]
    aux_fasta = Renome_seq(key,fastah)

    db = open('genome_pm_db.fa','w') 
    for acc in aux_fasta:
        db.write('>' + acc + '\n')
        db.write(str(aux_fasta[acc]) + '\n')
    db.close()
    
    query = open('query.fa','w')
    query.write('>'+key+'\n')
    query.write(str(seq) + '\n')
    query.close()

    #cmd = 'cat ' + sys.argv[2] +' >> genome_pm_db.fa'
    #process = subprocess.Popen(cmd,stdout=subprocess.PIPE, shell=True)
    #proc_stdout = process.communicate()[0].strip()

    cmd = 'lastdb -P96 -uNEAR -cR11 genome_db ' + 'genome_pm_db.fa\n'
    cmd += 'lastal -m100 -D10000000 -P100 -u1 genome_db query.fa | last-split > genome_aln.maf\n'
    cmd += 'maf-convert tab genome_aln.maf | awk -F\'=\' \'$2 <= 10e-7\' > ' + key.replace('|arrow','') + '.tab\n'
    print (cmd)
    comm = open ('command.bash','w')
    comm.write(cmd)
    comm.close()
    cmd = 'bash command.bash'

    process = subprocess.Popen(cmd,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    
    tabfile = key.replace('|arrow','') + '.tab'
    tabfile = tabfile.replace('.tab','_f10.tab')
    
    f = open(key.replace('|arrow','') + '.tab','r')
    f2 = open(tabfile,'w')
    
    count_hits = defaultdict(int)
    hits = defaultdict(list)
    for line in f:
        if '#' not in line:
            line = line.strip()
            l = line.split('\t')
            hits[l[1]].append(l)
            count_hits[l[1]] +=1
        else:
            f2.write(str(line))
    
    for k in sorted(count_hits):
        if int(count_hits[k]) > int(sys.argv[3]):
            size = len(hits.get(k))
            i=0
            while i < size:
                for t in hits[k][i]:
                    if t == '-':
                        f2.write ('-\t')
                    else:
                        f2.write (str(t)+'\t')
                i +=1
                f2.write('\n')
    f2.close()
    
    print (tabfile)
    f2 = open (tabfile,'r')  
    p_homologous = defaultdict(list)
    for line in f2:
        if '#' in line:
            continue
        else:
            l = line.split('\t')
            p_homologous[l[6]].append(l[1])

    ref_file = open('ref.fa','w')
    ref_file.write('>'+key+'\n')
    ref_file.write(str(seq) + '\n')
    ref_file.close()    
    
    query_file = open('qry.fa','w')
    
    for q in p_homologous:
        acc = list(dict.fromkeys(p_homologous[q]))
        for a in acc:
            if a in fasta:
            	query_file.write('>'+ a + '\n')
            	query_file.write(str(fasta[a]) + '\n')
            else:
            	query_file.write('>'+ a + '\n')
            	query_file.write(str(fastah[a]) + '\n')
    query_file.close()
    
    cmd = '/home/watson/Downloads/mummer-4.0.0beta2/nucmer -t 100 --maxmatch ref.fa -c 65 -l 500 qry.fa --prefix' + ' ' + key + '\n'
    cmd += 'delta-filter -i 0.7 -q ' + key.replace('|arrow','') + '.delta > ' + key.replace('|arrow','') +'.q.delta' + '\n'
    cmd += 'mummerplot --postscript ' + key.replace('|arrow','') +'.q.delta -p '+ key.replace('|arrow','') + '.q.delta -R ref.fa -Q qry.fa --filter --layout \n'
    print ("---------------------------------------" + key.replace('|arrow','') + "-----------------------------------------------")
    cmd_file = open('cmd.sh','w')
    cmd_file.write(cmd)
    cmd_file.close()

    cmd = 'bash cmd.sh'
    process = subprocess.Popen(cmd,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
