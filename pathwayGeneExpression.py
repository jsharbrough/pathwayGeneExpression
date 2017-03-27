def pathwayGOTerms(pathways,goTerms):
    infile = open(pathways,'r')
    infile2 = open(goTerms,'r')
    pathwayTranscriptDict = {}
    transcriptGOTermDict = {}
    transcriptDict = {}
    pathwayList = []
    for line in infile:
        line2 = line
        while line2[-1] == '\n' or line2[-1] == '\t' or line2[-1] == '\r':
            line2 = line2[0:-1]
        lineSplit = line2.split('\t')
        pathway = lineSplit[0]
        transcript = lineSplit[1]
        profile = int(lineSplit[2])
        if pathway in pathwayTranscriptDict:
            transcripts = pathwayTranscriptDict[pathway]
            transcripts.append(transcript)
            pathwayTranscriptDict[pathway] = transcripts
        else:
            pathwayTranscriptDict[pathway] = [transcript]
            pathwayList.append(pathway)
        if transcript not in transcriptDict:
            transcriptDict[transcript] = profile
    infile.close() 
    for line in infile2:
        line2 = line
        while line2[-1] == '\n' or line2[-1] == '\t' or line2[-1] == '\r':
            line2 = line2[0:-1]
        lineSplit2 = line2.split('\t')
        transcript2 = lineSplit2[1]
        goTerms = lineSplit2[2]
        goTermSplit = goTerms.split(';')
        transcriptGOTermDict[transcript2] = goTermSplit 
    infile2.close()
    GOPathwayDict = {}
    pathwayGODict = {}
    maxProfDict = {}
    pathwaysWithDominantProfile = []
    for pathway in pathwayList:
        transcripts = pathwayTranscriptDict[pathway]
        num1 = 0
        num2 = 0
        num3 = 0
        num4 = 0
        num5 = 0
        num6 = 0
        pathwayGoTerms = []
        for transcript in transcripts:
            expProfile = transcriptDict[transcript]
            if expProfile == 1:
                num1 += 1
            elif expProfile == 2:
                num2 += 1
            elif expProfile == 3:
                num3 += 1
            elif expProfile == 4:
                num4 += 1
            elif expProfile == 5:
                num5 += 1
            elif expProfile == 6:
                num6 += 1
            goTerms = transcriptGOTermDict[transcript]
            for term in goTerms:
                if term not in pathwayGoTerms:
                    pathwayGoTerms.append(term)
        pathwayGODict[pathway] = pathwayGoTerms
        profileNums = [num1,num2,num3,num4,num5,num6]
        maxProfiles = max(profileNums)
        listy = []
        i = 1
        for prof in profileNums:
            if prof == maxProfiles:
                listy.append(i)
            i += 1
        if len(listy) == 1 and maxProfiles > 1:
            maxProfDict[pathway] = listy[0]
            pathwaysWithDominantProfile.append(pathway)
            for goTerm in pathwayGoTerms:
                if goTerm in GOPathwayDict:
                    pathways = GOPathwayDict[goTerm]
                    if pathway not in pathways:
                        pathways.append(pathway)
                    GOPathwayDict[goTerm] = pathways
                else:
                    GOPathwayDict[goTerm] = [pathway]
    outfile = open('outfile.txt','w')
    outfile.write('GO_term\tPathways\t#_Observed_Shared_Expression_Profiles_Within_GO_Term\t#_Expected_Shared_Expression_Profiles_Within_GO_Term\t#_Possible_Shared_Expression_Profiles_Within_GO_Term\t%_Shared_Expression_Profiles_Within_GO_Term\t#_Observed_Shared_Expression_Profiles_Outside_GO_Term\t#_Expected_Shared_Expression_Profiles_Outside_GO_Term\t#_Possible_Shared_Expression_Profiles_Outside_GO_Term\t%_Shared_Expression_Profiles_Outside_GO_Term\n')
    for goTerm in GOPathwayDict:
        pathways = GOPathwayDict[goTerm]
        if len(pathways) > 1:
            numPossible = float(((len(pathways))*(len(pathways)-1))/2)
            numExpected = float(numPossible/6)
            numShared = 0
            for i in range(len(pathways)-1):
                profile1 = maxProfDict[pathways[i]]
                for j in range(i+1,len(pathways)):
                    profile2 = maxProfDict[pathways[j]]
                    if profile1 == profile2:
                        numShared += 1
            percentShared = round(float(numShared)/numPossible,2)
            numBetween = len(pathwaysWithDominantProfile) - len(pathways)
            numBetweenPossible = float((numBetween*(numBetween-1))/2)
            numBetweenExpected = float(numBetweenPossible/6)
            numSharedBetween = 0
            for pathway in pathways:
                profile1 = maxProfDict[pathway]
                for outsidePathway in pathwaysWithDominantProfile:
                    profile2 = maxProfDict[outsidePathway]
                    if outsidePathway not in pathways:
                        if profile1 == profile2:
                            numSharedBetween += 1
            percentSharedBetween = round(float(numSharedBetween)/numBetweenPossible,2)
            outfile.write(goTerm + '\t' + str(pathways) + '\t' + str(numShared) + '\t' + str(numExpected) + '\t' + str(numPossible) + '\t' + str(percentShared) + '\t' + str(numSharedBetween) + '\t' + str(numBetweenExpected) + '\t' + str(numBetweenPossible) + '\t' + str(percentSharedBetween) + '\n')
    outfile2 = open('pathwayProfiles.txt','w')
    print maxProfDict
    for pathway in pathwayList:
        outfile2.write(pathway + '\t')
        GoTerms = pathwayGODict[pathway]
        for term in GoTerms:
            outfile.write(term + '\t')
        if pathway in maxProfDict:
            outfile2.write(str(maxProfDict[pathway]))
        else:
            outfile2.write('N/A')
        for term in GoTerms:
            outfile2.write('\t' + term)
        outfile2.write('\n') 
    outfile.close()
    infile2.close()
    outfile2.close()
       
def withinPathwayVsBetweenPathwayExpression():
    infile = open('pathway_gene_homeolog_order.txt','r')
    outfile = open('Within_vs_Between_Shared_Expression_Exact_Matches_Genes.txt','w')
    outfile2 = open('Within_vs_Between_Shared_Expression_Exact_Matches_Pathways.txt','w')
    outfile.write('Gene' + '\t' + 'Pathway' + '\t' + 'Differentially_Expressed_Within_Pathway_Not_Shared_Profiles' + '\t' + 'Differentially_Expressed_Within_Pathway_Shared_Profiles' + '\t' + 'Differentially_Expressed_Between_Pathway_Not_Shared_Profiles' + '\t' + 'Differentially_Expressed_Between_Pathway_Shared_Profiles' + '\n')
    outfile2.write('Pathway' + '\t' + 'Number_Genes_In_Pathway' + '\t' + 'Differentially_Expressed_Within_Pathway_Not_Shared_Profiles' + '\t' + 'Differentially_Expressed_Within_Pathway_Shared_Profiles' + '\t' + 'Differentially_Expressed_Between_Pathway_Not_Shared_Profiles' + '\t' + 'Differentially_Expressed_Between_Pathway_Shared_Profiles' + '\n')
    geneList = []
    geneDict = {}
    pathwayList = []
    for line in infile:
        lineSplit = line.split('\t')
        pathway = lineSplit[0]
        gene = lineSplit[1]
        expression = lineSplit[2]
        geneList.append(gene)
        if pathway not in pathwayList:
            pathwayList.append(pathway)
        geneDict[gene] = (pathway,expression[0:-1])
    pathwayDict = {}
    for pathway in pathwayList:
        pathwayDict[pathway] = []
    for gene in geneList:
        geneInfo = geneDict[gene]
        pathway = geneInfo[0]
        expression = geneInfo[1]
        pathwayDict[pathway] += [(gene,expression)]
    for i in range(len(geneList)):
        gene1 = geneList[i]
        if i == 0:
            withinTotal = 0
            withinShared = 0
            betweenTotal = 0
            betweenShared = 0
            list1 = range(1,len(geneList))
            for j in list1:
                gene2 = geneList[j]
                gene1Info = geneDict[gene1]
                gene2Info = geneDict[gene2]
                pathway1  = gene1Info[0]
                pathway2  = gene2Info[0]
                expression1 = int(gene1Info[1])
                expression2 = int(gene2Info[1])
                if pathway1 == pathway2:
                    if expression1 == expression2:
                        withinTotal += 1
                        withinShared += 1
                    else:
                        withinTotal += 1
                else:
                    if expression1 == expression2:
                        betweenTotal += 1
                        betweenShared += 1
                    else:
                        betweenTotal += 1 
        elif i > 0  and i < (len(geneList)-1):
            withinTotal = 0
            withinShared = 0
            betweenTotal = 0
            betweenShared = 0
            list2 = range(0,i) + range(i+1,len(geneList))
            for j in list2:
                gene2 = geneList[j]
                gene1Info = geneDict[gene1]
                gene2Info = geneDict[gene2]
                pathway1  = gene1Info[0]
                pathway2  = gene2Info[0]
                expression1 = int(gene1Info[1])
                expression2 = int(gene2Info[1])
                if pathway1 == pathway2:
                    if expression1 == expression2:
                        withinTotal += 1
                        withinShared += 1
                    else:
                        withinTotal += 1
                else:
                    if expression1 == expression2:
                        betweenTotal += 1
                        betweenShared += 1
                    else:
                        betweenTotal += 1
        else:
            withinTotal = 0
            withinShared = 0
            betweenTotal = 0
            betweenShared = 0
            list3 = range(len(geneList)-1)
            for j in list3:
                gene2 = geneList[j]
                gene1Info = geneDict[gene1]
                gene2Info = geneDict[gene2]
                pathway1  = gene1Info[0]
                pathway2  = gene2Info[0]
                expression1 = int(gene1Info[1])
                expression2 = int(gene2Info[1])
                if pathway1 == pathway2:
                    if expression1 == expression2:
                        withinTotal += 1
                        withinShared += 1
                    else:
                        withinTotal += 1
                else:
                    if expression1 == expression2:
                        betweenTotal += 1
                        betweenShared += 1
                    else:
                        betweenTotal += 1
        if withinTotal > 0:
            percentSharedWithin = float(withinShared)/float(withinTotal)
            percentSharedBetween = float(betweenShared)/float(betweenTotal)
            outfile.write(gene1 + '\t' + pathway1 + '\t' + str(withinTotal - withinShared) + '\t' + str(withinShared) + '\t' + str(betweenTotal - betweenShared) + '\t' + str(betweenShared) + '\n')
    for path in pathwayList:
        pathInfo = pathwayDict[path]
        pathExpression = []
        genesInPath = []
        for gene in pathInfo:
            pathExpression.append(gene[1])
            genesInPath.append(gene[0])
        numPossibleWithinPath = len(pathInfo)*(len(pathInfo)-1)/2
        numPossibleBetweenPath = (len(geneList)-len(pathInfo))*(len(pathInfo))
        withinPathShared = 0
        betweenPathShared = 0
        if len(pathInfo) > 1:
            for i in range(len(pathInfo)-1):
                expValue1 = pathExpression[i]
                list1 = pathExpression[i+1:]
                for expValue2 in list1:
                    if expValue1 == expValue2:
                        withinPathShared += 1
                for gene in geneList:
                    outsideGeneInfo = geneDict[gene]
                    if outsideGeneInfo[0] != path:
                        expValue3 = outsideGeneInfo[1]
                        if expValue1 == expValue3:
                            betweenPathShared += 1
            percentSharedWithin = float(withinPathShared)/float(numPossibleWithinPath)
            percentSharedBetween = float(betweenPathShared)/float(numPossibleBetweenPath) 
            outfile2.write(path + '\t' + str(len(pathInfo)) + '\t' +str(numPossibleWithinPath - withinPathShared) + '\t' + str(withinPathShared) + '\t' + str(numPossibleBetweenPath - betweenPathShared) + '\t' + str(betweenPathShared) + '\n')       
    infile.close()
    outfile.close()
    outfile2.close()
        
    
def pathwayGeneExpression():
    infile = open('pathway_list.txt','r')
    infile1 = open('pathway_gene_homeolog_order.txt','r')
    pathwayList = []
    for line in infile:
        pathwayName = line[0:-1]
        if pathwayName not in pathwayList:
            pathwayList.append(pathwayName)
    infile.close()
    geneInfoList = []
    for line in infile1:
        lineSplit = line[0:-1].split('\t')
        pathway = lineSplit[0]
        gene = lineSplit[1]
        expressionOrder = lineSplit[2]
        geneInfoList.append([pathway,gene,expressionOrder[0]])
    pathwayMatches = {}
    for pathway in pathwayList:
        currList = []
        pathwayDict = {}
        for gene in geneInfoList:
            if gene[0] == pathway:
                pathwayDict[gene[1]] = gene[2]
                currList.append(gene[1])
        matches = 0
        highestOrMatches = 0
        lowestOrMatches = 0
        middleOrMatches = 0
        for i in range(len(currList)):
            geneA = currList[i]
            currExpressionOrder = pathwayDict[geneA]
            for j in range(1,len(currList)):
                geneB = currList[j]
                geneBExpressionOrder = pathwayDict[geneB]
                if geneA != geneB:
                    if currExpressionOrder == geneBExpressionOrder:
                        matches += 1
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif currExpressionOrder == '1':
                        if geneBExpressionOrder == '2':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '3':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '6':
                            middleOrMatches += 1
                    elif currExpressionOrder == '2':
                        if geneBExpressionOrder == '1':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '5':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '4':
                            middleOrMatches += 1
                    elif currExpressionOrder == '3':
                        if geneBExpressionOrder == '4':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '1':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '5':
                            middleOrMatches += 1
                    elif currExpressionOrder == '4':
                        if geneBExpressionOrder == '3':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '6':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '4':
                            middleOrMatches += 1
                    elif currExpressionOrder == '5':
                        if geneBExpressionOrder == '6':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '2':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '3':
                            middleOrMatches += 1
                    elif currExpressionOrder == '6':
                        if geneBExpressionOrder == '5':
                            highestOrMatches += 1
                        elif geneBExpressionOrder == '4':
                            lowestOrMatches += 1
                        elif geneBExpressionOrder == '1':
                            middleOrMatches += 1
                         
        numPossibleMatches = ((len(currList)*(len(currList)-1)))/2.0
        if len(currList) > 1:
            pathwayMatches[pathway] = [len(currList),numPossibleMatches/6.0,1/6.0, matches,numPossibleMatches, matches/numPossibleMatches,numPossibleMatches/3.0, 1/3.0, highestOrMatches, highestOrMatches/numPossibleMatches, lowestOrMatches, lowestOrMatches/numPossibleMatches, middleOrMatches, middleOrMatches/numPossibleMatches]
        else:
            pathwayMatches[pathway] = [len(currList),0,0,0,0,0,0,0,0,0,0,0,0,0]
    outfile = open('pathwayMatches.txt','w')
    outfile.write('Pathway' + '\t' + '# Diff Expressed Genes' + '\t' + '# Expected Matches' + '\t' + 'Expected Fraction Shared Expression' + '\t' + '# Matching Expression Orders' + '\t' + '# Possible Matches' + '\t' + 'Fraction Shared Expression' + '\t' + '# Expected Or Matches' + '\t' + '# Highest Genome Matches' + '\t' + 'Fraction Highest Genome Matches' + '\t' + '# Lowest Genome Matches' + '\t' + 'Fraction Lowest Genome Matches' + '\t' + '# Middle Genome Matches' + '\t' + 'Fraction Middle Genome Matches' + '\n')
    for pathway in pathwayList:
        currPath = pathwayMatches[pathway]
        outfile.write(pathway + '\t' + str(currPath[0]) + '\t' + str(currPath[1]) + '\t' + str(currPath[2]) + '\t' + str(currPath[3]) + '\t' + str(currPath[4]) + '\t' + str(currPath[5]) + '\t' + str(currPath[6]) + '\t' + str(currPath[7]) + '\t' + str(currPath[8]) + '\t' + str(currPath[9]) + '\t' + str(currPath[10]) + '\t' + str(currPath[11]) + '\t' + str(currPath[12]) + '\t' + str(currPath[13]) + '\n')
    infile1.close()
    outfile.close()

def transcriptomeExpressionOrder():
    infile1 = open('pathway_gene_homeolog_order.txt','r')
    geneList = []
    geneDict = {}
    for line in infile1:
        lineSplit = line[0:-2].split('\t')
        gene = lineSplit[1]
        geneList.append(gene)
        geneDict[gene] = lineSplit[2]
    matches = 0
    highestOrMatches = 0
    lowestOrMatches = 0
    middleOrMatches = 0
    for i in range(len(geneList)):
        geneA = geneList[i]
        currExpressionOrder = geneDict[geneA]
        for j in range(1,len(geneList)):
            geneB = geneList[j]
            geneBExpressionOrder = geneDict[geneB]
            if geneA != geneB:
                if currExpressionOrder == geneBExpressionOrder:
                    matches += 1
                elif currExpressionOrder == '1':
                    if geneBExpressionOrder == '1':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '2':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '3':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '6':
                        middleOrMatches += 1
                elif currExpressionOrder == '2':
                    if geneBExpressionOrder == '2':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '1':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '5':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '4':
                        middleOrMatches += 1
                elif currExpressionOrder == '3':
                    if geneBExpressionOrder == '3':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '4':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '1':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '5':
                        middleOrMatches += 1
                elif currExpressionOrder == '4':
                    if geneBExpressionOrder == '4':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '3':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '6':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '4':
                        middleOrMatches += 1
                elif currExpressionOrder == '5':
                    if geneBExpressionOrder == '5':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '6':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '2':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '3':
                        middleOrMatches += 1
                elif currExpressionOrder == '6':
                    if geneBExpressionOrder == '6':
                        highestOrMatches += 1
                        lowestOrMatches += 1
                        middleOrMatches += 1
                    elif geneBExpressionOrder == '5':
                        highestOrMatches += 1
                    elif geneBExpressionOrder == '4':
                        lowestOrMatches += 1
                    elif geneBExpressionOrder == '1':
                        middleOrMatches += 1
    numPossibleMatches = ((len(geneList)*(len(geneList)-1)))/2.0
    matchDetails = [len(geneList),numPossibleMatches/6.0,1/6.0, matches,numPossibleMatches, matches/numPossibleMatches,numPossibleMatches/3.0, 1/3.0, highestOrMatches, highestOrMatches/numPossibleMatches, lowestOrMatches, lowestOrMatches/numPossibleMatches, middleOrMatches, middleOrMatches/numPossibleMatches]
    outfile = open('pathwayMatches_total.txt','w')
    outfile.write('# Diff Expressed Genes' + '\t' + '# Expected Matches' + '\t' + 'Expected Fraction Shared Expression' + '\t' + '# Matching Expression Orders' + '\t' + '# Possible Matches' + '\t' + 'Fraction Shared Expression' + '\t' + '# Expected Or Matches' + '\t' + '# Highest Genome Matches' + '\t' + 'Fraction Highest Genome Matches' + '\t' + '# Lowest Genome Matches' + '\t' + 'Fraction Lowest Genome Matches' + '\t' + '# Middle Genome Matches' + '\t' + 'Fraction Middle Genome Matches' + '\n')
    outfile.write(str(matchDetails[0]) + '\t' + str(matchDetails[1]) + '\t' + str(matchDetails[2]) + '\t' + str(matchDetails[3]) + '\t' + str(matchDetails[4]) + '\t' + str(matchDetails[5]) + '\t' + str(matchDetails[6]) + '\t' + str(matchDetails[7]) + '\t' + str(matchDetails[8]) + '\t' + str(matchDetails[9]) + '\t' + str(matchDetails[10]) + '\t' + str(matchDetails[11]) + '\t' + str(matchDetails[12]) + '\t' + str(matchDetails[13]) + '\n')
    infile1.close()
    outfile.close()           
        