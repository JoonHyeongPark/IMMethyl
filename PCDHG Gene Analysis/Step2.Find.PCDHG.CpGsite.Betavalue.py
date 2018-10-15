target = "PCDHG"
file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice/"
input_target_site = open(file_path + target + ".CpGsites.Information.txt", 'r')

data = input_target_site.read().splitlines()
gene_index = data[0].split().index("UCSC_RefGene_Name")
target_site = list(map(lambda x : data[x].split()[0], range(1, len(data))))
target_gene = list(map(lambda x : data[x].split()[gene_index], range(1, len(data))))

cancerlist = ["PANCANCER"]

index = 0
for i in range(len(cancerlist)) :

    input_tumor = open(file_path + cancerlist[i] + ".humanmethylation450.tumor.31tumors.txt", 'r')
    input_normal = open(file_path + cancerlist[i] + ".humanmethylation450.normal.31tumors.txt", 'r')

    output_tumor = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.tumor.txt", 'w')
    output_normal = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.normal.txt", 'w')

    sample_tumor = input_tumor.readline().split()[1:]
    sample_normal = input_normal.readline().split()[1:]

    output_tumor.write("Site\tRelatedGene\t" + "\t".join(sample_tumor) + "\n")
    output_normal.write("Site\tRelatedGene\t" + "\t".join(sample_normal) + "\n")

    while(True) :

        line1 = input_tumor.readline().split()
        line2 = input_normal.readline().split()

        if(len(line2) == 0) : break

        if(line1[0] == target_site[index]) :

            index += 1
            output_tumor.write(line1[0] + "\t" + target_gene[index - 1] + "\t" + "\t".join(line1[1:]) + "\n")
            output_normal.write(line2[0] + "\t" + target_gene[index - 1] + "\t" + "\t".join(line2[1:]) + "\n")

            if(index == len(target_site)) : break
