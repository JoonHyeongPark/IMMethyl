# -*- coding: utf-8 -*-
cancerlist = ["PANCANCER", "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "UCEC", "UCS", "UVM"]
#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

########################################### 주의사항 : 완성된 파일에서 반드시 마지막 탭을 제거할 것 ##################################################

probe_count = 485577

######################################################################################################################################################

input_file1 = []
input_file2 = []

output_tumor1 = open("Excluded_Tumor_Cor_CpGsite_CytAct_pearson.txt", 'w')
output_tumor2 = open("Excluded_Tumor_Cor_CpGsite_CytAct_spearman.txt", 'w')

#tumor_header = "PANCANCER\tACC\tBLCA\tBRCA\tCESC\tCHOL\tCOAD\tDLBC\tESCA\tESCASTAD\tGBM\tHNSC\tKICH\tKIRC\tKIRP\tLGG\tLIHC\tLUAD\tLUSC\tMESO\tOV\tPAAD\tPCPG\tPRAD\tREAD\tSARC\tSKCM\tSTAD\tTGCT\tTHCA\tTHYM\tUCEC\tUCS\tUVM\n"
tumor_header = "\tPANCANCER_cor\tPANCANCER_P-value\tACC_cor\tACC_P-value\tBLCA_cor\tBLCA_P-value\tBRCA_cor\tBRCA_P-value\tCESC_cor\tCESC_P-value\tCHOL_cor\tCHOL_P-value\tCOAD_cor\tCOAD_P-value\tESCA_cor\tESCA_P-value\tESCASTAD_cor\tESCASTAD_P-value\tGBM_cor\tGBM_P-value\tHNSC_cor\tHNSC_P-value\tKICH_cor\tKICH_P-value\tKIRC_cor\tKIRC_P-value\tKIRP_cor\tKIRP_P-value\tLIHC_cor\tLIHC_P-value\tLUAD_cor\tLUAD_P-value\tLUSC_cor\tLUSC_P-value\tMESO_cor\tMESO_P-value\tPAAD_cor\tPAAD_P-value\tPCPG_cor\tPCPG_P-value\tPRAD_cor\tPRAD_P-value\tREAD_cor\tREAD_P-value\tSARC_cor\tSARC_P-value\tSKCM_cor\tSKCM_P-value\tSTAD_cor\tSTAD_P-value\tTGCT_cor\tTGCT_P-value\tTHCA_cor\tTHCA_P-value\tUCEC_cor\tUCEC_P-value\tUCS_cor\tUCS_P-value\tUVM_cor\tUVM_P-value\n"

for i in range(0, len(cancerlist)) :
    input_file1.append(open(cancerlist[i] + ".Tumor_Cor_CpGsite_CytAct_pearson.txt", 'r'))
    input_file2.append(open(cancerlist[i] + ".Tumor_Cor_CpGsite_CytAct_spearman.txt", 'r'))

    input_file1[i].readline()
    input_file2[i].readline()

output_tumor1.write(tumor_header); output_tumor2.write(tumor_header)

for site_number in range(0, probe_count) :
    
    value_print1 = ""; value_print2 = ""
    probe_name = "None"
    
    for i in range(0, len(cancerlist)) : 
        tumor1 = input_file1[i].readline().split()
        tumor2 = input_file2[i].readline().split()
        
        probe_name = tumor1.pop(0)
        del tumor2[0]
        
        for j in range(0, len(tumor1)) :
            value_print1 += "%s\t" % tumor1[j]
            value_print2 += "%s\t" % tumor2[j]
    
    output_tumor1.write("%s\t%s\n" % (probe_name, value_print1))
    output_tumor2.write("%s\t%s\n" % (probe_name, value_print2))

    if((site_number + 1) % 10000 == 0) : print(site_number + 1, " completed.\n")

for i in range(0, len(cancerlist)) :
    input_file1[i].close()
    input_file2[i].close()

output_tumor1.close()
output_tumor2.close()
        
print("END")
