import operator
import numpy

probe_count = 485577
#probe_count = 5

cancerlist = ["PANCANCER"]

input_file = []
output_file = []

sample_id = []

for k in range(0, len(cancerlist)) :

    input_file.append(open(cancerlist[k] + ".optimal_cutoff_data.txt", 'r'))
    output_file.append(open(cancerlist[k] + ".optimal_cutoff_result.txt", 'w'))

    sample_id.append(input_file[k].readline().split())

    for i in range(0, len(sample_id[k]) / 2) : del sample_id[k][i + 1]

sample_cutoff_count = [[0 for j in range(len(sample_id[i]))] for i in range(len(cancerlist))]

for k in range(0, len(cancerlist)) :

    for i in range(0, probe_count) :
    
        line = input_file[k].readline().split()

        if(len(line) == 0) :
            print(probe_count)
            continue
    
        del line[0]
        
        beta_value = []
    
        whole_tumor = 0; whole_normal = 0
    
        valid_data = 0

        for j in range(0, len(line) / 2) :
    
            if(line[j * 2] == "NA") : continue
    
            check = True
            if(line[j * 2 + 1] == "0") :
                check = False
                whole_normal += 1
    
            beta_value.append([float(line[j * 2]), check, j])
            valid_data += 1
    
        if(valid_data == 0) :
            continue
    
        whole_tumor = valid_data - whole_normal
    
        beta_value = sorted(beta_value, key = operator.itemgetter(0))
    
        count_tumor = []; count_normal = []
    
        count_tumor.append(beta_value[0][1])
        count_normal.append(1 - beta_value[0][1])
    
        for j in range(1, valid_data) :
            count_tumor.append(count_tumor[j - 1] + beta_value[j][1])
            count_normal.append(count_normal[j - 1] + (1 - beta_value[j][1]))
    
        cutoff = []
    
        tumor_low = []; tumor_high = []
        normal_low = []; normal_high = []
    
        fpr = []; tpr = []

        argmax_parameter = []

        for j in range(0, valid_data - 1) : 
    
            cutoff.append(float(beta_value[j][0] + beta_value[j + 1][0]) / 2)
    
            tumor_low.append(float(count_tumor[j])) # normal -> tumor -> false negative
            tumor_high.append(float(whole_tumor - count_tumor[j])) # tumor -> tumor -> true positive
    
            normal_low.append(float(count_normal[j])) # normal -> normal -> true negative
            normal_high.append(float(whole_normal - count_normal[j])) # tumor -> normal -> false positive
    
            fpr.append(normal_high[j] / (normal_high[j] + normal_low[j]))
            tpr.append(tumor_high[j] / (tumor_high[j] + tumor_low[j]))
    
            argmax_parameter.append(tpr[j] - fpr[j])
    
        optimal_index = numpy.argmax(argmax_parameter)
        optimal_threshold = cutoff[optimal_index]
        optimal_tpr = tumor_high[optimal_index]

        for j in range(optimal_index + 1, valid_data) :
            if(beta_value[j][1]) : sample_cutoff_count[k][beta_value[j][2]] += 1

        #printline = "%f %d" % (optimal_threshold, optimal_tpr)

        #print(printline)

    print(cancerlist[k] + " completed.")

for i in range(0, len(cancerlist)) :

    for j in range(0, len(sample_id[i])) :

        printline = sample_id[i][j]
        printline += "\t%d\n" % sample_cutoff_count[i][j]

        output_file[i].write(printline)

for k in range(0, len(cancerlist)) :
    input_file[k].close()
    output_file[k].close()
