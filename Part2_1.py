genomes = [[] ] * 5
dict = ["Reston_genome", "Sudan_genome", "TaiForest_genome", "Zaire_genome", "Bundibugyo_genome"]
dict2 = ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"]
for i in range(5):
    address = "BioProjectFiles/" + dict[i] + ".fasta"
    f = open(address, 'r')
    genome = f.read()
    genomes[i] = genome
    genomes[i] = genomes[i].split("\n",1)[1]
marburg_genes = [[]] * 7
f = open("BioProjectFiles/Marburg_Genes.fasta",'r')
marburgGenes = f.read()
indexes = []
for i in range(len(marburgGenes)):
    if marburgGenes[i] == ">":
        indexes.append(i)
for i in range(6):
    marburg_genes[i] = marburgGenes[indexes[i]: indexes[i+1] - 1]
    marburg_genes[i] = marburg_genes[i].split("\n",1)[1]




def aligner(ref, read):
    n = len(ref)
    m = len(read)
    match = 1
    mismatch = -2
    gap_open = -3
    dp = {}
    for i in range(n+1):
        dp[0,i] = 0
    for i in range(1,m+1):
        dp[i,0] = dp[i-1,0] + gap_open
    for i in range(1, n+1):
        print(i)
        for j in range(1,m+1):
            if ref[i-1] == read[j-1]:
                dp[j,i] = max(dp[j-1,i-1] + match, dp[j-1,i] + gap_open, dp[j,i-1] + gap_open)
            if ref[i-1] != read[j-1]:
                dp[j,i] = max(dp[j-1,i-1] + mismatch, dp[j-1,i] + gap_open, dp[j,i-1] + gap_open)
    maximum = 0
    endIndexi = 0
    for i in range(n+1):
        if dp[m,i] > maximum:
            maximum = dp[m,i]
            endIndexi = i

    indexi = endIndexi
    indexj = m
    while indexj > 0:
        if dp[indexj, indexi] == dp[indexj - 1, indexi] + gap_open:
            indexj -= 1
        elif dp[indexj, indexi] == dp[indexj, indexi - 1] + gap_open:
            indexi -= 1
        elif dp[indexj, indexi] == dp[indexj - 1, indexi - 1] + mismatch or dp[indexj, indexi] == dp[indexj - 1, indexi - 1] + match:
            indexi -= 1
            indexj -= 1
    startIndexi = indexi
    return startIndexi, endIndexi - 1
alignes = [[]]
print((aligner(genomes[3], marburg_genes[0])))