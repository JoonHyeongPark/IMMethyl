import scipy.stats as stat

target = "PCDHG"
file_path = "/home/kwangsookim_lab/joonhyeong_park/Practice/"

cancerlist = ["PANCANCER"]

for i in range(len(cancerlist)) :

    input_tumor = open(file_path + cancerlist[i] + "." + target + ".Score.Matching.Of.Each.Sample.tumor.txt", 'r')
    input_normal = open(file_path + cancerlist[i] + "." + target + ".Score.Matching.Of.Each.Sample.Normal.txt", 'r')

    input_tumor.readline()
    input_normal.readline()

    tumor_betavalue_score = []
    tumor_expression_score = []

    while(True) :

        line = input_tumor.readline().split()
        if(len(line) == 0) : break

        tumor_betavalue_score.append(float(line[1]))
        tumor_expression_score.append(float(line[2]))

    normal_betavalue_score = []
    normal_expression_score = []

    while(True) :

        line = input_normal.readline().split()
        if(len(line) == 0) : break

        normal_betavalue_score.append(float(line[1]))
        normal_expression_score.append(float(line[2]))

    print(stat.spearmanr(tumor_betavalue_score, tumor_expression_score))
    print(stat.spearmanr(normal_betavalue_score, normal_expression_score))
