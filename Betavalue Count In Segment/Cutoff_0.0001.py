#input_file = open("Betavalue_Array_Sorted.csv", 'r')
input_file = open("Betavalue_Array_CountSort.csv", 'r')

threshold_array = []
count_array = []

last_index = []

def MakeThresholdArray() :
        
    global threshold_array, count_array, last_index
            
    gap = [0.0001]
                
    for i in range(0, len(gap)) : 
                        
        cutoff = 0.0
        threshold_array.append([])
        count_array.append([])
        last_index.append(0)

        count = 0
        
        while(cutoff <= 1.0) :    
            threshold_array[i].append(cutoff)
            count_array[i].append(0)
            count += 1
            cutoff = gap[i] * float(count)
                                        
    return

def Process() :

    line = input_file.readline()
    debug_count = 1
            
    while(len(line) > 1) :
 
        value = float(line)
#        value = float(line.replace('"', '')[:7])

        for j in range(0, len(threshold_array)) : 
            for k in range(last_index[j], len(threshold_array[j]) - 1) :
                if(value < threshold_array[j][k + 1]) :
                    count_array[j][k] += 1
                    last_index[j] = k
                    break

        line = input_file.readline()
        debug_count += 1
                    
        if(debug_count % 10000 == 0) : print(str(debug_count) + " completed.")
                                                           
    for j in range(0, len(threshold_array)) :
                                                                   
        output_file = open(str(threshold_array[j][1]) + "cutoff.count_array.txt", 'w')
                                                                           
        for k in range(0, len(threshold_array[j]) - 1) :
            printline = str(threshold_array[j][k]) + " : " + str(count_array[j][k]) + "\n"
            output_file.write(printline)

    return

MakeThresholdArray()
Process()
