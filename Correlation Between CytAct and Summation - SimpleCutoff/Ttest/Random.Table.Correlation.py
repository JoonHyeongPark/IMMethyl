import random

cancerlist = ["PANCANCER"]

p_threshold = [0.05, 0.005, 0.0005, 0.0000001]
percentage = [0.01, 0.025, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]

CpGsites_number = []

test_type = "TTest"

sample_number = 8357
execution_number = 10000

for i in range(len(cancerlist)) :

    for j in range(len(p_threshold)) :

        input_cpgsites = open(str(p_threshold[j]) + "." + cancerlist[i] + ".CpGsites.By." + test_type + ".txt", 'r')
        input_cpgsites.readline() # junk line

        lines1 = input_cpgsites.read().splitlines()
        CpGsites_number.append(len(lines1))
       
        for k in range(len(percentage)) :

            input_samples = open("Pvalue." + str(p_threshold[j]) + ".Percentage." + str(percentage[k]) + "." + cancerlist[i] + "." + "Positive.MeaningfulCpGsites.By.Ttest_Binarization.Summation.By.Tstat.txt", 'r')
            lines2 = input_samples.read().splitlines()

            sample_id = []
            total_check = 0

            for line in lines2 : 
                line = line.split()
                sample_id.append(line[0])
                total_check += int(line[1])

            total = CpGsites_number[j] * sample_number

            print("%s %s, #CpGsites : %s, #Samples : %s, #Total : %s" % (str(p_threshold[j]), str(percentage[k]), str(CpGsites_number[j]), str(len(sample_id)), str(total_check)))
            print(total)

            for iteration in range(execution_number) :

                random_summation = []
                total_random = 0

                for l in range(sample_number) :

                    random_temp = random.randrange(0, CpGsites_number[j] + 1)
    
                    if(total_random + random_temp > total) : random_temp = random.randrange(0, total - total_random + 1)
                    elif(total_random + random_temp == total) : random_temp = 0
    
                    total_random += random_temp
                    random_summation.append(random_temp)
    
                random.shuffle(random_summation)

#                output.write(random_summation)
