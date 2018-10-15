import numpy

cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "STES", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "LAML"]

input_data_path = "/ngs/data3/public_data/TCGA_FireBrowse/mRNASeq"
merge_data_path = "/home/kwangsookim_lab/joonhyeong_park/Practice"

input_target = open(merge_data_path + "/Targeted.Genes.txt", 'r')
target = input_target.read().splitlines()

input_mRNA = []

output_merge_mRNA = open(merge_data_path + "/PANCANCER.PCDHG.Gene.mRNAseq_RSEM_normalized.log2.32tumors.txt", "w")

mRNA_header = ["Gene"]

sample_id = []
flat_id = []

for i in range(len(cancerlist)) :

    mRNA_path = input_data_path + "/" + cancerlist[i] + "/gdac.broadinstitute.org_" + cancerlist[i] + ".mRNAseq_Preprocess.Level_3.2016012800.0.0/" + cancerlist[i] + ".uncv2.mRNAseq_RSEM_normalized_log2.txt"
    input_mRNA.append(open(mRNA_path, 'r'))

    sample_id.append(input_mRNA[i].readline().replace('\n', '').split()[1:])
    flat_id += sample_id[i]

target_dic = {}
for each_id in flat_id : target_dic[each_id] = numpy.zeros((len(target)), dtype = float)
output_merge_mRNA.write("Gene\t" + "\t".join(flat_id) + "\n")

while(True) :

    empty_line_check = 0
    
    for i in range(len(cancerlist)) :
        
        line = input_mRNA[i].readline().split()
        
        length = len(line)
        empty_line_check += length

        if(length == 0) : continue

        line[0] = line[0][:line[0].index("|")]

        if(line[0] in target) :

            target_index = target.index(line[0])
            del line[0]

            for j in range(len(line)) :
                if(line[j] == "NA") : target_dic[sample_id[i][j]][target_index] = -1
                else : target_dic[sample_id[i][j]][target_index] = line[j]

    if(empty_line_check == 0) : break

A_TOTAL = {}
B_TOTAL = {}

for each_id in flat_id :
    A_TOTAL[each_id] = 0.0
    B_TOTAL[each_id] = 0.0

for i in range(len(target)) :
    output_merge_mRNA.write(target[i])

    for each_id in flat_id :
        if(target_dic[each_id][i] <= 0) :
            output_merge_mRNA.write("\tNA")
            continue
        output_merge_mRNA.write("\t" + str(target_dic[each_id][i]))

        if("PCDHGA" in target[i]) : A_TOTAL[each_id] += target_dic[each_id][i]
        else : B_TOTAL[each_id] += target_dic[each_id][i]

    output_merge_mRNA.write("\n")

output_merge_mRNA.write("PCDHGA_TOTAL")
for each_id in flat_id : output_merge_mRNA.write("\t%f" % A_TOTAL[each_id])

output_merge_mRNA.write("\nPCDHGA_MEAN")
for each_id in flat_id : output_merge_mRNA.write("\t%f" % (A_TOTAL[each_id] / 12))

output_merge_mRNA.write("\nPCDHGB_TOTAL")
for each_id in flat_id : output_merge_mRNA.write("\t%F" % B_TOTAL[each_id])

output_merge_mRNA.write("\nPCDHGB_MEAN")
for each_id in flat_id : output_merge_mRNA.write("\t%f" % (B_TOTAL[each_id] / 8))

output_merge_mRNA.write("\nPCDHG_SCORE")
for each_id in flat_id : output_merge_mRNA.write("\t%f" % ((A_TOTAL[each_id] + B_TOTAL[each_id]) / 20))
