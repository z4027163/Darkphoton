file1 = open('2017step2.log', 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
    line = line.strip('\n')
    line = line.strip('\t')
    count += 1
    x=line.split("/")
    y=x[1].split("_")
    print("    \"tree_%s   %s\""%(y[-1],x[0]))
    
