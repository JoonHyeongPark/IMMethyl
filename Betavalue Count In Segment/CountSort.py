input_file = open("Betavalue_Array.csv", 'r')
output_file = open("Betavalue_Array_CountSort.csv", 'w')

count_array = []

def Process() :

    line = input_file.readline()
    debug_count = 1
           
    for i in range(0, 10 ** 7 + 1) :
        count_array.append(0)

    while(len(line) > 1) :

        value = int(float(line) * (10 ** 7))
        count_array[value] += 1

        line = input_file.readline()
        debug_count += 1
                    
        if(debug_count % 10000 == 0) :
            print(str(debug_count) + " completed.")

    for i in range(0, 10 ** 7 + 1) :

        for j in range(0, count_array[i]) :
            printline = str(float(i) / float(10 ** 7)) + "\n"
            output_file.write(printline)
    return

Process()
