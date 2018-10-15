cancerlist = ["PANCANCER"]

probe_site = 485577

for i in range(len(cancerlist)) :

    input_file = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')
    output_file = open(cancerlist[i] + ".NA.sites.txt", 'w')

    input_file.readline()
    input_file.readline()

    for j in range(probe_site) :

        line = input_file.readline().split()
        site_id = line.pop(0)

        check = False
        for k in range(len(line)) :
            
            if(line[k] != "NA") :
                check = True
                break

        if(check == False) : output_file.write(site_id + "\n")
        if(j % 10000 == 0) : print("%d completed." % j)
