target = "PCDHG"
file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice/"

cancerlist = ["PANCANCER"]

for i in range(len(cancerlist)) :
    input_tumor = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.tumor.txt", 'r')
    input_normal = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.normal.txt", 'r')

    input_expression = open(file_path + cancerlist[i] + "." + target + ".Gene.mRNAseq_RSEM_normalized.log2.32tumors.txt", "r")

    output_tumor_expression = open(file_path + cancerlist[i] + "." + target + ".Expression.Pattern.tumor.txt", 'w')
    output_normal_expression = open(file_path + cancerlist[i] + "." + target + ".Expression.Pattern.normal.txt", 'w')

    sample_tumor = input_tumor.readline().split()[2:]
    sample_normal = input_normal.readline().split()[2:]

    sample_tumor = list(map(lambda x : sample_tumor[x][:15].replace('-', ''), range(len(sample_tumor))))
    sample_normal = list(map(lambda x : sample_normal[x][:15].replace('-', ''), range(len(sample_normal))))

    sample_expression = input_expression.readline().split()[1:]

    tumor_expression_index = []
    normal_expression_index = []

    tumor_header = "Gene"
    normal_header = "Gene"

    for j in range(len(sample_expression)) :
        check_sample = sample_expression[j][:15].replace('-', '')

        if(check_sample in sample_tumor) :
            tumor_expression_index.append(sample_tumor.index(check_sample))
            tumor_header += "\t" + sample_expression[j]

        else : tumor_expression_index.append(-1)

        if(check_sample in sample_normal) :
            normal_expression_index.append(sample_normal.index(check_sample))
            normal_header += "\t" + sample_expression[j]

        else : normal_expression_index.append(-1)

    output_tumor_expression.write(tumor_header + "\n")
    output_normal_expression.write(normal_header + "\n")

    while(True) :

        expression_line = input_expression.readline().split()
        if(len(expression_line) == 0) : break

        gene_name = expression_line.pop(0)

        output_tumor_expression.write(gene_name)
        output_normal_expression.write(gene_name)

        for j in range(len(sample_expression)) :

            if(tumor_expression_index[j] != -1) : output_tumor_expression.write("\t" + expression_line[j])
            if(normal_expression_index[j] != -1) : output_normal_expression.write("\t" + expression_line[j])

        output_tumor_expression.write("\n")
        output_normal_expression.write("\n")
