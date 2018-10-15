division_part = 20
execution_number = 500

cancerlist = ["PANCANCER"]

for k in range(len(cancerlist)) :

    output_file = open("RANDOM.MERGE." + cancerlist[k] + ".Correlation.Values." + str(execution_number * division_part) + ".Times.txt", 'w')

    for i in range(division_part) :

        number = ""

        if(i < 9) : number = "0" + str(i + 1)
        else : number = str(i + 1)

        input_file = open("RANDOM." + number + "." + cancerlist[k] + ".Correlation.Values." + str(execution_number) + ".Times.txt", 'r')

        lines = input_file.read().splitlines()
        for line in lines : output_file.write(line + "\n")
