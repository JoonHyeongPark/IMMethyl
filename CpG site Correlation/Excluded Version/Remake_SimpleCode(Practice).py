# -*- coding: utf-8 -*-
from scipy import stats

cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PRAD"]

probe_count = 485577
#probe_count = 10

pan_betavalue = []
pan_cytact = []
sample_id = []

sample_table = {}

input_file1 = []
output_file_pearson = []
output_file_spearman = []

input_cytact = open("TCGA_methylation_cowork_1.txt",'r')

header = input_cytact.readline().split()
indexing = header.index("CytAct")

whole_table = input_cytact.readlines()

for line in whole_table :
    line = line.split()
    ID = line[0].replace("_", "")

    sample_table[ID] = line[indexing]

sample_name = []

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
    output_file_pearson.append(open(cancerlist[i] + ".pearson.txt", 'w'))
    output_file_spearman.append(open(cancerlist[i] + ".spearman.txt", 'w'))

    sample_name.append(input_file1[i].readline().split())

    del sample_name[i][0]; del sample_name[i][0]

    for j in range(0, len(sample_name[i])) : sample_name[i][j] = sample_name[i][j][:15].replace('-', "")

    input_file1[i].readline()

output_whole_pearson = open("whole_pearson.txt", 'w')
output_whole_spearman = open("whole_spearman.txt", 'w')

output_pan_pearson = open("PAN.pearson.txt", 'w')
output_pan_spearman = open("PAN.spearman.txt", 'w')

for probe_iteration in range(0, probe_count) :

    pearson_printline = ""
    spearman_printline = ""

    probe_name = ""

    for i in range(0, len(cancerlist)) :

        line = input_file1[i].readline().split()

        probe_name = line.pop(0)
    
        betavalue_vector = []
        cytact_vector = []

        for k in range(0, len(line)) :  
            if(sample_name[i][k] in sample_table and line[k] != "NA") : 
                betavalue_vector.append(float(line[k]))
                cytact_vector.append(float(sample_table[sample_name[i][k]]))
                
                pan_betavalue.append(float(line[k]))
                pan_cytact.append(float(sample_table[sample_name[i][k]]))

        cor1 = stats.pearsonr(betavalue_vector, cytact_vector)
        cor2 = stats.spearmanr(betavalue_vector, cytact_vector)

        each_pearson = "%s\t%s\t" % (str(cor1[0]), str(cor1[1]))
        each_spearman = "%s\t%s\t" % (str(cor2[0]), str(cor2[1]))

        pearson_printline += "%s" % each_pearson
        spearman_printline += "%s" % each_spearman

        output_file_pearson[i].write(probe_name); output_file_pearson[i].write("\t"); output_file_pearson[i].write(each_pearson); output_file_pearson[i].write("\n")
        output_file_spearman[i].write(probe_name); output_file_spearman[i].write("\t"); output_file_spearman[i].write(each_spearman); output_file_spearman[i].write("\n")

    cor1 = stats.pearsonr(pan_betavalue, pan_cytact)
    cor2 = stats.spearmanr(pan_betavalue, pan_cytact)

    pan_pearson = probe_name + "\t%s\t%s\t" % (str(cor1[0]), str(cor1[1]))
    pan_spearman = probe_name + "\t%s\t%s\t" % (str(cor2[0]), str(cor2[1]))

    pearson_printline = pan_pearson + pearson_printline
    spearman_printlnie = pan_spearman + spearman_printline

    output_whole_pearson.write(pearson_printline); output_whole_pearson.write("\n")
    output_whole_spearman.write(spearman_printline); output_whole_spearman.write("\n")

    output_pan_pearson.write(pan_pearson); output_pan_pearson.write("\n")
    output_pan_spearman.write(pan_spearman); output_pan_spearman.write("\n")

    if(probe_iteration % 100 == 0) : print(probe_iteration)
