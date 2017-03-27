def salamanderGenes():
    infile = open('salamander_BLAST.txt','r')
    allHitsList = []
    geneList = []
    geneDict = {}
    threeGeneList = []
    geneGenomeList = []
    for line in infile:
        if line[0] != '#':
            tabbedLine = line.split('\t')
            gene = tabbedLine[2]
            gene = gene.strip('\n')
            allHitsList.append(gene)
            genome = tabbedLine[1]
            geneGenome = (gene,genome)
            geneGenomeList.append(geneGenome)
            if gene not in geneList:    
                geneList.append(gene)
    for gene in geneList:
        geneDict[gene] = 0
    for gene in geneList:
        geneCounter = 0
        for item in allHitsList:
            if gene == item:
                geneCounter += 1
        geneDict[gene] = geneCounter
    for gene in geneDict:
        if geneDict[gene] == 3:
            threeGeneList.append(gene)
    genomeDict = {}
    for gene in threeGeneList:
        genomeDict[gene] = []
    for item in geneGenomeList:
        for gene in threeGeneList:
            if gene == item[0]:
                a = genomeDict[gene]
                a.append(item[1])
                genomeDict[gene] = a
    genesOfInterest = []
    for locus in genomeDict:
        if 'AL1' in genomeDict[locus]:
            if 'A1237' in genomeDict[locus]:
                if 'A1238' in genomeDict[locus]:
                    genesOfInterest.append(locus)
    outfile_AL1 = open('salamander_genes_AL1.fas','w')
    outfile_A1237 = open('salamander_genes_A1237.fas','w')
    outfile_A1238 = open('salamander_genes_A1238.fas','w')
    for gene in genesOfInterest:
        outfile_AL1.write('>' + gene + '\n')
        outfile_A1237.write('>' + gene + '\n')
        outfile_A1238.write('>' + gene + '\n')
    outfile_AL1.close()
    outfile_A1237.close()
    outfile_A1238.close()
    infile.close()

def salamanderPullSequences():
    infile_1 = open('AL1_annotated.fasta','r')
    infile_2 = open('salamander_genes_AL1.fas','w')
    
            