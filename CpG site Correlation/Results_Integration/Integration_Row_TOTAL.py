# -*- coding: utf-8 -*-
target_cancer = "PANCANCER"

input_file1 = []
input_file2 = []

output_tumor1 = open(target_cancer + "_Tumor_Cor_CpGsite_CytAct_pearson.txt", 'w')
output_tumor2 = open(target_cancer + "_Tumor_Cor_CpGSite_CytAct_spearman.txt", 'w')

for i in range(0, 8) :

    name = str(int(i + 1))

    input_file1.append(open(target_cancer + ".SEP_" + name + ".Tumor_Cor_CpGsite&CytAct_pearson.txt", 'r'))
    input_file2.append(open(target_cancer + ".SEP_" + name + "..Tumor_Cor_CpGSite&CytAct_spearman.txt", 'r'))

    if(i == 0) :
        output_tumor1.write(input_file1[i].readline())
        output_tumor2.write(input_file2[i].readline())

    else :
        input_file1[i].readline()
        input_file2[i].readline()

    printarray1 = input_file1[i].readlines()
    printarray2 = input_file2[i].readlines()

    for line1 in printarray1 : output_tumor1.write(line1)
    for line2 in printarray2 : output_tumor2.write(line2)

output_tumor1.close()
output_tumor2.close()
        
print("END")
