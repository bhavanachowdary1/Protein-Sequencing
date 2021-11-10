"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f=open(filename)
    text=f.read()
    textsplit=text.splitlines()
    string = "".join(textsplit) 
    return string
   


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    lstofcodons=[]
    replacedlst=[]
    for i in range(startIndex,len(dna),3):
        lstofcodons.append(dna[i:i+3])
        if dna[i:i+3]=='TAG' or dna[i:i+3]=='TAA' or dna[i:i+3]=='TGA':
            break
    for string in lstofcodons:
        string=string.replace("T","U")
        replacedlst.append(string)
    return replacedlst

'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f=open(filename)
    file=json.load(f)
    dictionary={}
    for i,j in file.items():
        for k in j:
            a=k.replace("T","U")
            dictionary[a]=i
    return dictionary



'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
    
def generateProtein(codons, codonD):
    protein=[]
    if codons[0]=='AUG':
        protein.append('Start')
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            protein.append(codonD[codons[i]])
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    a=readFile(dnaFilename)
    #print(a)
    b=makeCodonDictionary(codonFilename)
    i=0
    count=0
    proteinlist=[]
    while i<len(a):
        if a[i:i+3] == "ATG":
            c=dnaToRna(a,i)
            d=generateProtein(c,b)
            proteinlist.append(d)
            i=i+3*len(c)
        else:
            i=i+1
            count=count+1
    print("Unused bases are : ", count)
    print("Total Bases are : ", len(a)) 
    print("Proteins count are : ", len(proteinlist))
    return proteinlist


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    uniqueproteinlist=[]
    for i in proteinList1:
        for j in proteinList2:
            if i==j and i not in uniqueproteinlist:
                uniqueproteinlist.append(i)
    #print(uniqueproteinlist)
    return uniqueproteinlist



'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combinedlist=[]
    for i in proteinList:
        for j in i:
            combinedlist.append(j)
    return combinedlist


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dictionary={}
    for i in aaList:
        if i not in dictionary:
            dictionary[i]=1
        else:
            dictionary[i]+=1
    return dictionary


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    lst1=combineProteins(proteinList1)
    lst2=combineProteins(proteinList2)
    dictionary1=aminoAcidDictionary(lst1)
    dictionary2=aminoAcidDictionary(lst2)
    freq1={}
    freq2={}
    temp=[] #without appending start stop 
    final_difference=[]
    for i in dictionary1:
        freq1[i]=dictionary1[i]/len(lst1)
        # print(freq1[i])
        if i not in temp and i!='Start' and i!='Stop':
            temp.append(i)
    #print("this is i===============================",temp)
    for j in dictionary2:
        freq2[j]=dictionary2[j]/len(lst2)
        if j not in temp and j!='Start' and j!='Stop':
            temp.append(j)
    #print("this is j===============================",temp)
    for k in temp:
        frequency1=0
        frequency2=0
        if k in freq1:
            frequency1=freq1[k]
        if k in freq2:
            frequency2=freq2[k]
        diff=frequency2-frequency1
        if diff < -cutoff or diff > cutoff:
            difflst=[k,frequency1,frequency2]
            final_difference.append(difflst)
    #print(final_difference)
    return final_difference



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        cproteins=""
        lst=i[1:(len(i)-1)]
        count=0
        for j in lst:
            cproteins+=j
            count+=1
            if count!=len(lst):
                cproteins+="-"
        if len(cproteins)!=0:
            print(cproteins)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in differences:
        a=i[0]
        freq1=round(i[1]*100,2)
        freq2=round(i[2]*100,2)
        print(str(a)+" "+str(freq1)+" % in Seq1"+","+str(freq2)+"% in Seq2")
    return



def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    gene=[]
    lst1=combineProteins(proteinList1)
    lst2=combineProteins(proteinList2)
    dictionary1=aminoAcidDictionary(lst1)
    dictionary2=aminoAcidDictionary(lst2)
    for i in dictionary1:
        if i not in gene:
            gene.append(i)
    for j in dictionary2:
        if j not in gene:
            gene.append(j)
    gene.sort()
    return gene
     #aminoacids sorted list in gene



  


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    lst=combineProteins(proteinList)
    dictionary=aminoAcidDictionary(lst)
    newlist=[] # frequencies of the labels
    for i in labels:
        if i in dictionary:
            a=dictionary[i]/len(lst)
            newlist.append(a)
        else:
            newlist.append(0)
    #print(newlist)
    return newlist



'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    w=0.4
    xvalues=np.arange(len(xLabels))
    plt.bar(xvalues,freqList1,width=-w,align='edge',label=label1,edgecolor=edgeList)
    plt.bar(xvalues,freqList2,width=w,align='edge',label=label2,edgecolor=edgeList)
    plt.xticks(ticks=list(range(len(xLabels))),labels=xLabels,rotation="horizontal")
    plt.legend()
    plt.title("Comparision of Frequencies")
    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    edgelst=[]
    lst=[]
    for i in range(len(biggestDiffs)):
        lst.append(biggestDiffs[i][0]) #0 element in 0 index 3 element list
    for i in range(len(labels)):
        if labels[i] in lst:
            edgelst.append("black")
        else:
            edgelst.append("white")
    return edgelst



'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    import matplotlib.pyplot as plt
    humanproteins=synthesizeProteins("data/human_p53.txt","data/codon_table.json")
    elephantproteins=synthesizeProteins("data/elephant_p53.txt","data/codon_table.json")
    cproteins=commonProteins(humanproteins,elephantproteins)
    diff=findAminoAcidDifferences(humanproteins,elephantproteins,0.005)
    displayTextResults(cproteins,diff)
    labels=makeAminoAcidLabels(humanproteins,elephantproteins)
    f1=setupChartData(labels,humanproteins)
    f2=setupChartData(labels,elephantproteins)
    edges=makeEdgeList(labels,diff)
    createChart(labels, f1, "Human", f2, "Elephant", edgeList=edges)

    
    #print(labels)
    return





### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    #test.testDnaToRna()
    #test.testMakeCodonDictionary()
    #test.testGenerateProtein()
    #test.testReadFile()
    #test.testSynthesizeProteins()
    #test.testCommonProteins()
    #test.testCombineProteins()
    #test.testAminoAcidDictionary()
    #test.testFindAminoAcidDifferences()
    #test.testMakeAminoAcidLabels()
    #test.testSetupChartData()
    #test.testCreateChart()
    #test.testMakeEdgeList()



    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
