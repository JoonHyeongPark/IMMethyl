# -*- coding: utf-8 -*-
cancerlist = ["ESCASTAD", "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
#cancerlist = ["STAD", "OV"]
########################################### 주의사항 : 완성된 파일에서 반드시 마지막 탭을 제거할 것 ##################################################

probe_count = 485577
sample_count = 0

######################################################################################################################################################

input_file1 = []
#input_file2 = []

output_tumor = open("PANCANCER.humanmethylation450.tumor.txt", 'w')
#output_normal = open("PANCANCER.humanmethylation450.normal.txt", 'w')

tumor_header = "Hybridization\tREF\t"; #normal_header = "Hybridization\tREF\t"

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
#    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))
    
    sample_name1 = input_file1[i].readline().split()
#    sample_name2 = input_file2[i].readline().split()
    
    del sample_name1[0]; del sample_name1[0]
#    del sample_name2[0]; del sample_name2[0]
    
    for j in range(0, len(sample_name1)) : tumor_header += "%s\t" % sample_name1[j]    
#    for j in range(0, len(sample_name2)) : normal_header += "%s\t" % sample_name2[j]
        
    input_file1[i].readline(); #input_file2[i].readline() # 쓰레기 line 제거

output_tumor.write(tumor_header); #output_normal.write(normal_header)

output_tumor.write("\n############################################################################################################\n")
#output_normal.write("\n############################################################################################################\n")

for site_number in range(0, probe_count) :
    
    value_print1 = ""; #value_print2 = ""
    probe_name = "None"
    
    for i in range(0, len(cancerlist)) : 
        tumor = input_file1[i].readline().split()
        #normal = input_file2[i].readline().split()
        
        probe_name = tumor.pop(0); #probe_name = normal.pop(0)
        
        for j in range(0, len(tumor)) : value_print1 += "%s\t" % tumor[j]
        #for j in range(0, len(normal)) : value_print2 += "%s\t" % normal[j]
    
    output_tumor.write("%s\t%s\n" % (probe_name, value_print1)); #output_normal.write("%s\t%s\n" % (probe_name, value_print2))

    if((site_number + 1) % 10000 == 0) : print(site_number + 1, " completed.\n")

for i in range(0, len(cancerlist)) :
    input_file1[i].close()
#    input_file2[i].close()
    
output_tumor.close()
#output_normal.close()
        
print("END")
