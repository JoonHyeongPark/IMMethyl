cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_data_path = "/ngs/data3/public_data/TCGA_FireBrowse/Methylation"
merge_data_path = "/home/kwangsookim_lab/joonhyeong_park/Practice"

input_tumor = []
input_normal = []

output_merge_tumor = open(merge_data_path + "/PANCANCER.humanmethylation450.tumor.31tumors.txt", "w")
output_merge_normal = open(merge_data_path + "/PANCANCER.humanmethylation450.normal.31tumors.txt", 'w')

tumor_header = ["Site"]
normal_header = ["Site"]

for i in range(len(cancerlist)) :

    cancer_data_path = input_data_path + "/" + cancerlist[i] + ".humanmethylation450/" + cancerlist[i] + ".humanmethylation450."

    input_tumor.append(open(cancer_data_path + "tumor.txt", 'r'))
    input_normal.append(open(cancer_data_path + "normal.txt", 'r'))

    tumor_header += input_tumor[i].readline().replace('\n', '').split()[2:]
    input_tumor[i].readline()

    normal_header += input_normal[i].readline().replace('\n', '').split()[2:]
    input_normal[i].readline()

output_merge_tumor.write("\t".join(tumor_header) + "\n")
output_merge_normal.write("\t".join(normal_header) + "\n")

iteration_number = 0

while(True) :

    empty_line_check = 0
    for i in range(len(cancerlist)) :
        
        tumor_line = input_tumor[i].readline().split()
        normal_line = input_normal[i].readline().split()

        if(len(normal_line) == 0) : break

        empty_line_check += len(normal_line) - 1

        if(i != 0) :
            del tumor_line[0]
            del normal_line[0]

        output_merge_tumor.write("\t".join(tumor_line))
        output_merge_normal.write("\t".join(normal_line))

        if(i != len(cancerlist) - 1) :
            output_merge_tumor.write("\t")
            output_merge_normal.write("\t")

    if(empty_line_check == 0) : break

    output_merge_tumor.write("\n")
    output_merge_normal.write("\n")

    if(iteration_number % 10000 == 0) : print(iteration_number)
    iteration_number += 1
