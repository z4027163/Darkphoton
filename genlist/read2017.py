file1 = open('tem2017C.log', 'r')
Lines1 = file1.readlines()
for line in Lines1:
    line = line.strip('\n')
    line = line.strip('\t')
    x=line.split("/")
    y=x[1].split("_")
    print("    \"scout_%s  2017/ScoutingRunC/ScoutingCaloMuon/crab_20200617_140145/200617_120348/0000/ 2017C\""%(y[-1]))

file2 = open('tem2017D.log', 'r')
Lines2 = file2.readlines()
for line in Lines2:
    line = line.strip('\n')
    line = line.strip('\t')
    x=line.split("/")
    y=x[1].split("_")
    print("    \"scout_%s  2017/ScoutingRunD/ScoutingCaloMuon/crab_20200617_140351/200617_120401/0000/ 2017D\""%(y[-1]))

file3 = open('tem2017E.log', 'r')
Lines3 = file3.readlines()
for line in Lines3:
    line = line.strip('\n')
    line = line.strip('\t')
    x=line.split("/")
    y=x[1].split("_")
    print("    \"scout_%s  2017/ScoutingRunE/ScoutingCaloMuon/crab_20200617_140405/200617_120422/0000/ 2017E\""%(y[-1]))

file4 = open('tem2017F.log', 'r')
Lines4 = file4.readlines()
for line in Lines4:
    line = line.strip('\n')
    line = line.strip('\t')
    x=line.split("/")
    y=x[1].split("_")
    print("    \"scout_%s  2017/ScoutingRunF/ScoutingCaloMuon/crab_20200617_140425/200617_120436/0000/ 2017F\""%(y[-1]))
