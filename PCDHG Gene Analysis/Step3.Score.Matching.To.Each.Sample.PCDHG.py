import numpy

target = "PCDHG"
file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice/"

cancerlist = ["PANCANCER"]

for i in range(len(cancerlist)) :
    input_tumor = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.tumor.txt", 'r')
    input_normal = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.normal.txt", 'r')

    output_tumor = open(file_path + cancerlist[i] + "." + target + ".Score.Matching.Of.Each.Sample.tumor.txt", 'w')
    output_normal = open(file_path + cancerlist[i] + "." + target + ".Score.Matching.Of.Each.Sample.Normal.txt", 'w')

    header = "Sample\tBetavalueMean\tExpressionMean\n"
    output_tumor.write(header)
    output_normal.write(header)

    sample_tumor = input_tumor.readline().split()[2:]
    sample_normal = input_normal.readline().split()[2:]

    original_sample_tumor = []
    original_sample_normal = []

    tumor_sample_information = []
    normal_sample_information = []

    for j in range(len(sample_tumor)) :
        original_sample_tumor.append(sample_tumor[j])
        sample_tumor[j] = sample_tumor[j][:15].replace('-', '')
        tumor_sample_information.append({"Betavalue" : [], "Expression" : []})

    for j in range(len(sample_normal)) :
        original_sample_normal.append(sample_normal[j])
        sample_normal[j] = sample_normal[j][:15].replace('-','')
        normal_sample_information.append({"Betavalue" : [], "Expression" : []})

    input_expression = open(file_path + cancerlist[i] + "." + target + ".Gene.mRNAseq_RSEM_normalized.log2.32tumors.txt", 'r')
    sample_expression = input_expression.readline().split()[1:]

    tumor_expression_index = []
    normal_expression_index = []

    for j in range(len(sample_expression)) :

        check_sample = sample_expression[j][:15].replace('-', '')
        if(check_sample in sample_tumor) : tumor_expression_index.append(sample_tumor.index(check_sample))
        else : tumor_expression_index.append(-1)
        if(check_sample in sample_normal) : normal_expression_index.append(sample_normal.index(check_sample))
        else : normal_expression_index.append(-1)

    while(True) :

        tumor_line = input_tumor.readline().split()[2:]
        normal_line = input_normal.readline().split()[2:]

        if(len(normal_line) == 0) : break

        for j in range(len(tumor_line)) :
            if(tumor_line[j] == "NA") : continue
            tumor_sample_information[j]["Betavalue"].append(float(tumor_line[j]))

        for j in range(len(normal_line)) :
            if(normal_line[j] == "NA") : continue
            normal_sample_information[j]["Betavalue"].append(float(normal_line[j]))

    while(True) :

        expression_line = input_expression.readline().split()[1:]
        if(len(expression_line) == 0) : break

        for j in range(len(expression_line)) :

            if(tumor_expression_index[j] != -1) : 
                if(expression_line[j] == "NA") : tumor_sample_information[tumor_expression_index[j]]["Expression"].append(0)
                else : tumor_sample_information[tumor_expression_index[j]]["Expression"].append(float(expression_line[j]))

            if(normal_expression_index[j] != -1) :
                if(expression_line[j] == "NA") : normal_sample_information[normal_expression_index[j]]["Expression"].append(0)
                else : normal_sample_information[normal_expression_index[j]]["Expression"].append(float(expression_line[j]))

    for j in range(len(sample_tumor)) :

        if(len(tumor_sample_information[j]["Expression"]) != 0) :
            printline = original_sample_tumor[j] + "\t" + str(numpy.mean(tumor_sample_information[j]["Betavalue"])) + "\t" + str(tumor_sample_information[j]["Expression"][-1])
            output_tumor.write(printline + "\n")

    for j in range(len(sample_normal)) :

        if(len(normal_sample_information[j]["Expression"]) != 0) :
            printline = original_sample_normal[j] + "\t" + str(numpy.mean(normal_sample_information[j]["Betavalue"])) + "\t" + str(normal_sample_information[j]["Expression"][-1])
            output_normal.write(printline + "\n")
