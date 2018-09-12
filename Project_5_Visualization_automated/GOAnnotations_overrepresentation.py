import re
import os

#goAnnotationsReformatted="C:/Users/deluca/data/GO/go_annotations.txt"
annos = {}
cwd = os.getcwd()
os.chdir(cwd)
goAnnotations="mgi.gaf"
goOntology = "go.obo"

print("Loading GO Annotations")
file = open(goAnnotations, "r") 
for line in file:
    if line[0] != '!' :
        split = line.split('\t')
        gene = split[2]
        go = split[4]
        
        geneList = annos.setdefault(go,set())
        geneList.add(gene)
print(" .... loaded")

#with open(goAnnotationsReformatted,"w") as f:
#    for go, genes in annos.items():
#        f.write(go)
#        f.write('\t')
#        f.write(','.join(genes))
#        f.write('\n')
#

print("Loading GO Ontology")
ontology = {}
lastID=""
file = open(goOntology, "r") 
for line in file:
    if line.startswith("id:"):
        lastID= line.split(": ")[1][:-1]
        ontology.setdefault(lastID,[])
    elif line.startswith("is_a"):
        parent = line.split(" ")[1]
        ontology[lastID].append(parent)
print(" .... loaded")


print("Reversing GO Ontology")

ontologyR = {}
for child,parents in ontology.items():
    ontologyR.setdefault(child,[]) # the children need their own nodes, even if they will have no childen
    for parent in parents:
        node = ontologyR.setdefault(parent,[])
        node.append(child)

print(" .... reversed")


def getGenes(term,ontologyR,annos):
    genes = set()
    if term not in ontologyR:
        print(term,"not found in ontology")
        return genes
    if term in annos:
        genes = annos[term] 
    for child in ontologyR[term]:
        genes.update(getGenes(child,ontologyR,annos))
    return genes
        
        

#pantherResults = "/home/binliu/R/Project_3_Lung Fibrosis/Visualization_B_Tcells/Tcells_result/analysis_Tcells_overrepresentation_logFC<=-1_pvalue<=0.005.txt"
#pantherResultsAnno = "/home/binliu/R/Project_3_Lung Fibrosis/Visualization_B_Tcells/Tcells_result/analysis_Tcells_overrepresentation_logFC<=-1_pvalue<=0.005_Annotation.txt"

pantherResults = "Panther_result.txt"
pantherResultsAnno = os.path.join(cwd,"result/GO_to_annotation.txt")


print("Annotating panther results")
out = open(pantherResultsAnno,"w")
out.write("GO Term\tSet Size\tDE Genes\tExpected\tDirection\tEnrichment\tp-value\tAnnotations\tGenes\n")
file = open(pantherResults, "r")
i = 0
for line in file:
    i+=1
    if i > 13:
        split = line.split('\t')
        go = re.split("\(|\)",split[0])[1] # the first cell needs to be split into term and ID
        out.write(line[:-1])
        out.write('\t')
        # get the genes for this go term
        genes = getGenes(go,ontologyR,annos)
        out.write(str(len(genes)))
        out.write('\t')
        out.write(','.join(genes))
        out.write('\n')

file.close()
out.close()




