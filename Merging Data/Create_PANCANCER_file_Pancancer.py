#-*- coding: utf-8 -*-
cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

probe_count = 485577
sample_count = 0

input_file1 = []

output_tumor = open("PANCANCER.humanmethylation450.tumor.txt", 'w')

tumor_header = "Hybridization\tREF\t"; #normal_header = "Hybridization\tREF\t"

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
    
    sample_name1 = input_file1[i].readline().split()
    
    del sample_name1[0]; del sample_name1[0]
    
    for j in range(0, len(sample_name1)) : tumor_header += "%s\t" % sample_name1[j]    
        
    input_file1[i].readline(); #input_file2[i].readline() # 쓰레기 line 제거

output_tumor.write(tumor_header); #output_normal.write(normal_header)
output_tumor.write("\n############################################################################################################\n") # 형식 통일

for site_number in range(0, probe_count) :
    
    value_print1 = ""
    probe_name = "None"
    
    for i in range(0, len(cancerlist)) : 
        tumor = input_file1[i].readline().split()
        probe_name = tumor.pop(0)
        
        for j in range(0, len(tumor)) : value_print1 += "%s\t" % tumor[j]
    
    output_tumor.write("%s\t%s\n" % (probe_name, value_print1))

    if((site_number + 1) % 10000 == 0) : print(site_number + 1, " completed.\n")

for i in range(0, len(cancerlist)) :
    input_file1[i].close()
    
output_tumor.close()
        
print("END")
