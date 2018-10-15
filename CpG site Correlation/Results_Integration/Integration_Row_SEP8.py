# -*- coding: utf-8 -*-
target_cancer = "PANCANCER"

input_file1 = []
input_file2 = []

output_tumor1 = open(target_cancer + ".SEP_8.Tumor_Cor_CpGsite&CytAct_pearson.txt", 'w')
output_tumor2 = open(target_cancer + ".SEP_8..Tumor_Cor_CpGSite&CytAct_spearman.txt", 'w')

for i in range(0, 10) :

    name = str(i)

    input_file1.append(open(target_cancer + ".SEP_8." + name + ".Tumor_Cor_CpGsite&CytAct_pearson.txt", 'r'))
    input_file2.append(open(target_cancer + ".SEP_8." + name + "..Tumor_Cor_CpGSite&CytAct_spearman.txt", 'r'))

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

last1 = open(target_cancer + ".SEP_8._.Tumor_Cor_CpGsite&CytAct_pearson.txt", 'r')
last2 = open(target_cancer + ".SEP_8._..Tumor_Cor_CpGSite&CytAct_spearman.txt", 'r')

last1.readline(); last2.readline()

lastprint1 = last1.readlines(); lastprint2 = last2.readlines()

for line1 in lastprint1 : output_tumor1.write(line1)
for line2 in lastprint2 : output_tumor2.write(line2)

output_tumor1.close()
output_tumor2.close()
        
print("END")
