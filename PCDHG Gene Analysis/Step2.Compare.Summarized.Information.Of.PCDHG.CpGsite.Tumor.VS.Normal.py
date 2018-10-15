import numpy

target = "PCDHG"
file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice/"

cancerlist = ["PANCANCER"]
for i in range(len(cancerlist)) :

    input_tumor = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.tumor.txt", 'r')
    input_normal = open(file_path + cancerlist[i] + "." + target + ".Methylation.Pattern.normal.txt", 'r')

    header1 = input_tumor.readline().split()
    header2 = input_normal.readline().split()

    output = open(file_path + cancerlist[i] + "." + target + ".Compare.Methylation.Information.txt", 'w')
    output.write("Site\tRelatedGenes\tTumorMean\tNormalMean\tMeanGap\tTumorMedian\tNormalMedian\tMedianGap\tTumorMax\tNormalMax\tTumorMin\tNormalMin\n")

    while(True) : 

        line1 = input_tumor.readline().split()
        line2 = input_normal.readline().split()

        if(len(line2) == 0) : break

        raw_tumor_betavalue = line1[2:]
        raw_normal_betavalue = line2[2:]

        valid_tumor_betavalue = []
        valid_normal_betavalue = []

        for value in raw_tumor_betavalue :
            if(value == "NA") : continue
            valid_tumor_betavalue.append(float(value))

        for value in raw_normal_betavalue :
            if(value == "NA") : continue
            valid_normal_betavalue.append(float(value))

        output.write("\t".join(line1[:2]) + "\t")
        value = []

        if(len(valid_normal_betavalue) == 0 or len(valid_tumor_betavalue) == 0) :
            output.write("NotEnoughValidValues\n")
            continue

        value.append(str(numpy.mean(valid_tumor_betavalue)))
        value.append(str(numpy.mean(valid_normal_betavalue)))
        value.append(str(abs(float(value[0]) - float(value[1]))))

        value.append(str(numpy.median(valid_tumor_betavalue)))
        value.append(str(numpy.median(valid_normal_betavalue)))
        value.append(str(abs(float(value[3]) - float(value[4]))))

        value.append(str(max(valid_tumor_betavalue)))
        value.append(str(max(valid_normal_betavalue)))

        value.append(str(min(valid_tumor_betavalue)))
        value.append(str(min(valid_normal_betavalue)))

        output.write("\t".join(value) + "\n")
