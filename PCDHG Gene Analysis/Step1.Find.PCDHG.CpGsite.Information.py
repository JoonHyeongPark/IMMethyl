target = "PCDHG"

file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice"

input_file = open(file_path + "/HumanMethylation450_15017482_v1-2.prc.txt", 'r')
output_file = open(file_path + "/" + target + ".CpGsites.Information.txt", 'w')

header = input_file.readline().split()
gene_index = header.index("UCSC_RefGene_Name")

output_file.write("\t".join(header) + "\n")

while(True) :

    line = input_file.readline().split()

    if(len(line) == 0) : break

    gene_list = line[gene_index].split(';')

    check = False
    for gene in gene_list :

        if(target in gene) :
            check = True
            break

    if(check == True) : output_file.write("\t".join(line) + "\n")
