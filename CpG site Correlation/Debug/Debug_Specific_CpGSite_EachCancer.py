# -*- coding: utf-8 -*-
from scipy import stats

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PRAD"]
cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "UCEC", "UCS", "UVM", "PRAD"]

probe_count = 485577
target_site = "cg00066663"

betavalue_vector = []
cytact_vector = []
sample_id = []

sample_table = {}
check_table = {}

input_file1 = []

input_cytact = open("TCGA_methylation_cowork_1.txt",'r')

header = input_cytact.readline().split()
indexing = header.index("CytAct")

whole_table = input_cytact.readlines()

for line in whole_table :
    line = line.split()
    ID = line[0].replace("_", "")

    sample_table[ID] = line[indexing]

output = open("debug1.txt",'w')

output.write("%s\n" % target_site)

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
   
    sample_name = input_file1[i].readline().split()

    del sample_name[0]; del sample_name[0]

    for j in range(0, len(sample_name)) :
        sample_name[j] = sample_name[j][:15].replace("-", "") 

    input_file1[i].readline()

    for j in range(0, probe_count) :
        line = input_file1[i].readline().split()

        comparison = line.pop(0)

        if(comparison == target_site) :

            for k in range(0, len(line)) :  
                if(sample_name[k] in sample_table and line[k] != "NA" and sample_name[k] not in sample_id) : 
                    sample_id.append(sample_name[k])
                    betavalue_vector.append(float(line[k]))
                    cytact_vector.append(float(sample_table[sample_name[k]]))
            break

    input_file1[i].close()

for i in range(0, len(betavalue_vector)) :
    printline = "%s\t%s\t%s\n" % (sample_id[i], betavalue_vector[i], cytact_vector[i])
    output.write(printline)

cor = stats.spearmanr(betavalue_vector, cytact_vector)

lastline = "%f\t%f\n" % (cor[0], cor[1])

output.write(lastline)

output.close()
print("END")
