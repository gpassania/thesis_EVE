# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:04:04 2016

@author: zwhitfield
"""
import sys
import pandas as pd

outputdir = str(sys.argv[1])
fileName = str(sys.argv[2])

#Taxonomy categories to keep track of
taxonomyCategories = ['superkingdom','order','family','genus','species']

def LoadData (directory, fileNameWithEVEs, fileNameWithEVEhierarchies, fileNameWithEVEranks):

    #Read in potential files a bit differently given differing names/numbers of columns
    if 'closestTEtoEVEs' in fileNameWithEVEs:
        allEVEs = pd.read_csv(directory + fileNameWithEVEs,
                              names = ["ContigEVE","EVEstart","EVEend","EVEdescription","EVEscore","EVEstrand","EVEsomething","ContigTE","TEstart","TEend","TEdescription","TEscore","TEstrand","TEfamily","Distance"],
                              sep = "\t")
    elif 'Aag_Contigs' in fileNameWithEVEs:
        allEVEs = pd.read_csv(directory + fileNameWithEVEs,
                          names = ["ContigEVE","EVEstart","EVEend","EVEdescription","EVEscore","EVEstrand","EVEpi"],
                          sep = "\t")
    # elif 'TEsClosestToEVEs' in fileNameWithEVEs:
    #     allEVEs = pd.read_csv(directory + fileNameWithEVEs,
    #                       sep = "\t")

    EVEvirusHierarchies = pd.read_csv(directory + fileNameWithEVEhierarchies,
                              names = ["Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7", "Rank8", "Rank9"],
                              sep = "\t")
    
    EVERankHierarchies = pd.read_csv(directory + fileNameWithEVEranks,
                              names = ["A","B","C","D","E","F","G", "H", "I"],
                              sep = "\t")
    #All 3 variables above (allEVEs, EVEvirusHierarchies, EVERankHierarchies) should have the same number of lines. If not, something is wrong.

    #This organization is VERY specific to the NCBI virus taxonomy categories
    #This will be filled out as a dataframe filled with virus taxonomy assignments for each virus (one per line) in allEVEs.
    EVEhierarchiesOrganized = pd.DataFrame(index=range(0,len(EVEvirusHierarchies)), columns=taxonomyCategories)
    
    currentRow = 0
    for index, Rankrow in EVERankHierarchies.iterrows():#use pandas to iterate through each row of taxonomy ranks (again, each row corresponds to virus
        currentCol = 0
        for Rankcol in Rankrow: #iterate through each column
            #print EVErow[currentCol]
            if Rankcol in taxonomyCategories: #If current taxonomy rank is in taxonomyCategories, then use current row and Rank as coordinates as assign taxonomy
                EVEhierarchiesOrganized.ix[currentRow][Rankcol] = EVEvirusHierarchies.iloc[currentRow][currentCol] #This requires that EVEvirusHierarchies and EVERankHierarchies are exactly 1:1 (they should be)
            currentCol = currentCol + 1
        currentRow = currentRow + 1    

    #Merge files based on index, so now have taxonomy assignments for each virus in allEVEs
    Merged = pd.merge(allEVEs, EVEhierarchiesOrganized, right_index=True, left_index=True)
    
    #Hand annotate some viruses with missing families :(
    #Most of these are based on Li et al. eLife 2015,"Unprecedented genomic diversity of RNA viruses in arthropods reveals the ancestry of negative-sense RNA viruses"
    #Some are from Uniprot taxonomy categories
    #Not comprehensive changes; only those viruses I find in the EVE files
    Merged.loc[Merged.species == 'Wuhan Mosquito Virus 8', 'family'] = "Chuviridae"
    Merged.loc[Merged.species == 'Wuchang Cockraoch Virus 3', 'family'] = "Chuviridae"    
    Merged.loc[Merged.species == 'Lishi Spider Virus 1', 'family'] = "Chuviridae"    
    Merged.loc[Merged.species == 'Shayang Fly Virus 1', 'family'] = "Chuviridae"    
    Merged.loc[Merged.species == 'Wenzhou Crab Virus 2', 'family'] = "Chuviridae"    
    
    Merged.loc[Merged.species == 'Bole Tick Virus 2', 'family'] = "Rhabdoviridae"    
    Merged.loc[Merged.species == 'Shayang Fly Virus 2', 'family'] = "Rhabdoviridae"    
    Merged.loc[Merged.species == 'Wuhan Ant Virus', 'family'] = "Rhabdoviridae"    
    Merged.loc[Merged.species == 'Wuhan Fly Virus 2', 'family'] = "Rhabdoviridae"    
    Merged.loc[Merged.species == 'Wuhan House Fly Virus 1', 'family'] = "Rhabdoviridae"        
    Merged.loc[Merged.species == 'Wuhan Mosquito Virus 9', 'family'] = "Rhabdoviridae"
    Merged.loc[Merged.species == 'Yongjia Tick Virus 2', 'family'] = "Rhabdoviridae"
    
    Merged.loc[Merged.species == 'Cilv-C', 'family'] = "Virgaviridae"
    Merged.loc[Merged.species == 'Citrus leprosis virus C', 'family'] = "Virgaviridae"
    Merged.loc[Merged.species == 'Blueberry necrotic ring blotch virus', 'family'] = "Virgaviridae"
    
    Merged.loc[Merged.species == 'Wutai Mosquito Virus', 'family'] = "Bunyaviridae"

    Merged.to_csv(outputdir + fileNameWithEVEs + "_withTaxonomy.txt", sep='\t', header = True, index = False, quoting = False)
    
    return Merged
    
LoadData (outputdir,
          fileName,
          "TaxonomyHierarchyByEVE_" + fileName + ".txt",
          "RankingHierarchyByEVE_" + fileName + ".txt"
          )
