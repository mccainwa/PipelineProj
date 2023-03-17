import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#os.makedirs("PipelineProject_Walter_McCain") #Making New directory
os.chdir("PipelineProject_Walter_McCain") #Get into directory

#os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030'")
#os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033'")
#os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'")
#os.system("wget 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'")

#os.system("fastq-dump -I --split-files SRR5660030") #Uncompress downloaded data
#os.system("fastq-dump -I --split-files SRR5660033")
#os.system("fastq-dump -I --split-files SRR5660044")
#os.system("fastq-dump -I --split-files SRR5660045")

Entrez.email = 'wmccain@luc.edu'
handle = Entrez.efetch(db="nucleotide", id="NC_006273.2", rettype = "gb", retmode = "text")
record = SeqIO.read(handle, format = "genbank") #Read and seperate

outfile = open("fastafile.fasta", "a")

cds = []
for feature in record.features:
    if feature.type == 'CDS':
        cds.append(feature.extract(record.seq))
        recordid=feature.qualifiers["protein_id"]
        outfile.write(">" + str(recordid) + "\n")
        outfile.write(str(feature.extract(record.seq)) + "\n")

samples = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

for sample in samples:
    os.system(f"kallisto index -i {sample}.idx {sample}_1.fastq {sample}_2.fastq")
    os.system(f"kallisto quant -i {sample}.idx -o {sample}_output -t 4 {sample}_1.fastq {sample}_2.fastq")
    with open(f"{sample}_output/abundance.tsv", "w") as f:
        f.write("target_id\ttest_stat\tpval\tqval\n")

os.system("Rscript -e 'library(sleuth); metadata = data.frame(sample = c('SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045'), condition = c('2dpi', '2dpi', '6dpi', '6dpi')); so = sleuth_prep(metadata, ~condition, extra_bootstrap_summary=T); so = sleuth_fit(so); so = sleuth_fit(so, 'full'); res = sleuth_results(so, '2dpi', '6dpi'); write.table(res[res$qval<0.05,c('target_id','test_stat','pval','qval')], file='significant_transcripts.txt', sep='\t', row.names=F, quote=F)'")

with open("log.txt", "w") as f:
    f.write("target_id\ttest_stat\tpval\tqval\n")
with open("significant_transcripts.txt") as st:
    next(st)
    for line in st:
        fields = line.strip().split("\t")
        target_id = fields[0]
        test_stat = fields[1]
        pval = fields[2]
        qval = fields[3]
        f.write(f"{target_id}\t{test_stat}\t{pval}\t{qval}\n")

tpm_data = {}
for sample in samples:
    with open(f"{sample}_output/abundance.tsv") as f:
        f.readline() # skip header
        for line in f:
            fields = line.strip().split("\t")
            protein_id = fields[0]
            tpm = float(fields[4])
            if protein_id not in tpm_data:
                tpm_data[protein_id] = {}
            tpm_data[protein_id][sample] = tpm

output_data = []
for protein_id, sample_tpm in tpm_data.items():
    row = [protein_id]
    for sample in samples:
        row.append(sample_tpm.get(sample, 0))
    sample_tpm_values = list(sample_tpm.values())
    row.append(min(sample_tpm_values))
    row.append(sorted(sample_tpm_values)[len(sample_tpm_values)//2])
    row.append(sum(sample_tpm_values)/len(sample_tpm_values))
    row.append(max(sample_tpm_values))
    output_data.append(row)

with open("log.txt", "w") as f:
    f.write("protein_id\tsample1\tsample2\tsample3\tsample4\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
    for row in output_data:
        f.write("\t".join(str(x) for x in row) + "\n")
