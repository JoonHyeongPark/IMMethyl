input_file1 = open("Excluded_Tumor_Cor_CpGsite_CytAct_pearson.txt", 'r')
input_file2 = open("Excluded_Tumor_Cor_CpGsite_CytAct_spearman.txt", 'r')

probe_count = 485577
cancer_count = 33 # except PANCANCER

per = 0.8
percentage = float(cancer_count) * per

cutoff = 0.3

pearson_output1 = open(str(per) + "_pearson_" + str(cutoff) + "compact_version.txt", 'w')
pearson_output2 = open(str(per) + "_pearson_" + "-" + str(cutoff) + "compact_version.txt", 'w')
spearman_output1 = open(str(per) + "_spearman_" + str(cutoff) + "compact_versino.txt", 'w')
spearman_output2 = open(str(per) + "_spearman_" + "-" + str(cutoff) + "compact_version.txt", 'w')

header = input_file1.readline()
header = input_file2.readline()
pearson_output1.write(header)
pearson_output2.write(header)
spearman_output1.write(header)
spearman_output2.write(header)

for i in range(0, probe_count) :

    pearson_line = input_file1.readline().split()
    spearman_line = input_file2.readline().split()

    printline1 = pearson_line[0]
    printline2 = spearman_line[0]
    printline3 = pearson_line[0]
    printline4 = spearman_line[0]

    del pearson_line[0]; del spearman_line[0]

    p_pan_cor = pearson_line.pop(0)
    p_pan_P = pearson_line.pop(0)

    sp_pan_cor = spearman_line.pop(0)
    sp_pan_P = spearman_line.pop(0)

    if(p_pan_cor != "NA") :

        if(float(p_pan_P) < 0.05) :

            if(float(p_pan_cor) >= cutoff) : 

                count_3 = 0
                count_5 = 0
                count_7 = 0

                j = 0
                while j < cancer_count * 2 :

                        if(pearson_line[j] != "NA") :
                            if(float(pearson_line[j]) >= 0.7 and float(pearson_line[j+1]) < 0.05) : count_7 += 1
                            elif(float(pearson_line[j]) >= 0.5 and float(pearson_line[j+1]) < 0.05) : count_5 += 1
                            elif(float(pearson_line[j]) >= 0.3 and float(pearson_line[j+1]) < 0.05) : count_3 += 1

                        j += 2
    
                total_count = count_3 + count_5 + count_7
                if(total_count >= percentage) :
                    print("IN")
                    printline1 += "\t" + str(p_pan_cor) + "\t" + str(p_pan_P) + "\t"
                    for adding in pearson_line : printline1 += adding + "\t"
                    printline1 += "\n"

                    pearson_output1.write(printline1)


            if(float(p_pan_cor) <= 0 - cutoff) : 
    
                count_3 = 0
                count_5 = 0
                count_7 = 0
    
                j = 0
                while j < cancer_count * 2 :
    
                        if(pearson_line[j] != "NA") : 
                            if(float(pearson_line[j]) <= -0.7 and float(pearson_line[j+1]) < 0.05) : count_7 += 1
                            elif(float(pearson_line[j]) <= -0.5 and float(pearson_line[j+1]) < 0.05) : count_5 += 1
                            elif(float(pearson_line[j]) <= -0.3 and float(pearson_line[j+1]) < 0.05) : count_3 += 1
    
                        j += 2
    
                total_count = count_3 + count_5 + count_7
    
                if(total_count >= percentage) :
                    printline2 += "\t" + str(p_pan_cor) + "\t" + str(p_pan_P) + "\t"
                    for adding in pearson_line : printline2 += adding + "\t"
                    printline2 += "\n"

                    pearson_output2.write(printline2)


    if(sp_pan_cor != "NA") :

        if(float(sp_pan_P) < 0.05) :
    
            if(float(sp_pan_cor) >= cutoff) : 
    
                count_3 = 0
                count_5 = 0
                count_7 = 0

                j = 0
                while j < cancer_count * 2 :
  
                        if(spearman_line[j] != "NA") :

                            if(float(spearman_line[j]) >= 0.7 and float(spearman_line[j+1]) < 0.05) : count_7 += 1
                            elif(float(spearman_line[j]) >= 0.5 and float(spearman_line[j+1]) < 0.05) : count_5 += 1 
                            elif(float(spearman_line[j]) >= 0.3 and float(spearman_line[j+1]) < 0.05) : count_3 += 1
    
                        j += 2
    
                total_count = count_3 + count_5 + count_7

                if(total_count >= percentage) :
                    printline3 += "\t" + str(sp_pan_cor) + "\t" + str(sp_pan_P) + "\t"
                    print("IN")
                    for adding in spearman_line : printline3 += adding + "\t"
                    printline3 += "\n"

                    spearman_output1.write(printline3)

            if(float(sp_pan_cor) <= 0 - cutoff) : 
    
                count_3 = 0
                count_5 = 0
                count_7 = 0

                j = 0
                while j < cancer_count * 2 :
    
                        if(spearman_line[j] != "NA") : 
                            if(float(spearman_line[j]) <= -0.7 and float(spearman_line[j+1]) < 0.05) : count_7 += 1
                            elif(float(spearman_line[j]) <= -0.5 and float(spearman_line[j+1]) < 0.05) : count_5 += 1
                            elif(float(spearman_line[j]) <= -0.3 and float(spearman_line[j+1]) < 0.05) : count_3 += 1
    
                        j += 2
    
                total_count = count_3 + count_5 + count_7
    
                if(total_count >= percentage) : 
                    printline4 += "\t" + str(sp_pan_cor) + "\t" + str(sp_pan_P) + "\t"
                    print("IN")
                    for adding in spearman_line : printline4 += adding + "\t"
                    printline4 += "\n"

                    spearman_output2.write(printline4)
