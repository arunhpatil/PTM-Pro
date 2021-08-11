#!/usr/bin/python
reqList = list()
with open("ptms.txt", "r") as infile:
    inf = infile.readlines()
    for i in inf:
        if "PeptideSequence" not in i:
            i = i.strip()
            i = i.replace(" ","")
            i = i.split("\t")
            if "||" in i[3]:
                pipsitesplit = i[2].split("||")
                pipgenesplit = i[3].split("||")
                pipwinsplit = i[4].split("||")
                pipnovsplit = i[5].split("||")
                for index, pipe in enumerate(pipsitesplit):
#['YISATEGADsHPASGGGtPTGDEGQIGQVTTK', 'Phospho;Phospho', 'S162; T170||S162; T170', 'Sap47||Sap47', 'SATEGADsHPASGGG;HPASGGGtPTGDEGQ||SATEGADsHPASGGG;HPASGGGtPTGDEGQ', 'Yes;Yes||Yes;Yes']
                    if ";" in pipe:
                        pipemods = i[1].split(";")
                        pipesite = pipe.split(";")
                        pipewins = pipwinsplit[index].split(";")
                        pipenovt = pipnovsplit[index].split(";")
                        for scin, sc in enumerate(pipemods):
                            #print(i[0],sc,pipesite[scin],pipgenesplit[index],pipewins[scin],pipenovt[scin])
                            reqList.append([i[0],sc,pipesite[scin],pipgenesplit[index],pipewins[scin],pipenovt[scin]])
                        #semiKey = i[0]+"#"+i[1]
                    else:
                        reqList.append([i[0],i[1],pipe,pipgenesplit[index],pipwinsplit[index],pipnovsplit[index]])
                        #print(i[0],i[1],pipe,pipgenesplit[index],pipwinsplit[index],pipnovsplit[index])
                        #pass
            else:
                if ";" in i[1]:
                    semiptm = i[1].split(";")
                    semisite = i[2].split(";")
                    semiwin = i[4].split(";")
                    seminovelt = i[5].split(";")
                    for incx, site in enumerate(semiptm):
                        #print(i[0],site, semisite[incx],i[3],semiwin[incx],seminovelt[incx])
                        reqList.append([i[0],site, semisite[incx],i[3],semiwin[incx],seminovelt[incx]])
                        #pass
                else:
                    #print(i[0],i[1],i[2],i[3],i[4],i[5])
                    reqList.append([i[0],i[1],i[2],i[3],i[4],i[5]])

res = []
[res.append(x) for x in reqList if x not in res]
#print(len(reqList), len(res))
countList = list()
for i in res:
    site = i[2][0]
    #print(i[1], site, i[5])
    countList.append(str(i[1])+"_"+str(site)+"_"+str(i[5]))
    #countList.append([i[1], site, i[5]])
    #countList.append[str(i[1])+"_"+str(site)+str(i[5])]
    #countList.append[str(i[1])+"_"+str(site)+str(i[5])]
    #print(i[2][0])
    #print(i)
#for j in countList:
#    print(j)
duplicate_dict = {i:countList.count(i) for i in countList}
#print(duplicate_dict)
for k, y in duplicate_dict.items():
    print(k,y)
#['PeptideSequence', 'PTM_Modification', 'PTM_Site (Protein)', 'GeneSymbol', 'PTM_Window(PTM-Pro2.0)', 'dbPTM Evidence']
#['GASDsPAASR', 'Phospho', 'S45', 'rdgA', 'NRRGASDsPAASRAR', 'No']
#['RPVGGkPAAGGQR', 'Acetyl', 'K234', 'Ref1', 'FKRPVGGkPAAGGQR', 'No']
#['QLATkAAR', 'Acetyl', 'K24||K24', 'His3.3B||His3:CG33803', 'PRKQLATkAARKSAP||PRKQLATkAARKSAP', 'Yes||No']
#['GLGkGGAkR', 'Acetyl;Acetyl', 'K13; K17', 'His4:CG33909', 'KGGKGLGkGGAKRHR;GLGKGGAkRHRKVLR', 'Yes;Yes']

#['RLsAESPISEEDPVAATK', 'Phospho', 'S666||S629', 'Esyt2||Esyt2', 'QASSKRLsAESPISE||QASSKRLsAESPISE', 'Yes||Yes']
