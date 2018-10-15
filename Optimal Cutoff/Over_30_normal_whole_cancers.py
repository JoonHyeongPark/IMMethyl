# -*- coding: utf-8 -*-
#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PRAD"]
cancerlist = ["KIRC", "BRCA", "THCA", "HNSC", "LIHC", "PRAD", "UCEC", "KIRP", "LUSC", "COAD", "LUAD"]

probe_count = 485577 

input_file1 = []
input_file2 = []

tumor_sample_name = []
normal_sample_name = []

header = "\t"

output_file = open("PANCANCER.optimal_cutoff_data.txt", 'w')

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))

    tumor_sample_name.append(input_file1[i].readline().split())
    normal_sample_name.append(input_file2[i].readline().split())

    del tumor_sample_name[i][0]; del tumor_sample_name[i][0]
    del normal_sample_name[i][0]; del normal_sample_name[i][0]

    for j in range(0, len(tumor_sample_name[i])) :
        tumor_sample_name[i][j] = tumor_sample_name[i][j][:15].replace('-', "")
        header += "%s\tTumor/Normal\t" % tumor_sample_name[i][j]

    for j in range(0, len(normal_sample_name[i])) :
        normal_sample_name[i][j] = normal_sample_name[i][j][:15].replace('-', "")
        header += "%s\tTumor/Normal\t" % normal_sample_name[i][j]

    input_file1[i].readline()
    input_file2[i].readline()

header += "\n"
output_file.write(header)

for i in range(0, probe_count) :

    probe_name = ""
    printline = ""

    for j in range(0, len(cancerlist)) :

        line1 = input_file1[j].readline().split()
        line2 = input_file2[j].readline().split()

        probe_name = line1.pop(0)
        probe_name = line2.pop(0)

        for k in range(0, len(line1)) : printline += "%s\t%d\t" % (line1[k], j + 1)
        for k in range(0, len(line2)) : printline += "%s\t0\t" % line2[k]

    printline = probe_name + "\t" + printline + "\n"
    output_file.write(printline)

    if(i % 1000 == 0) : print(i)

output_file.close()

for i in range(0, len(cancerlist)) :
    input_file1[i].close()
    input_file2[i].close()
