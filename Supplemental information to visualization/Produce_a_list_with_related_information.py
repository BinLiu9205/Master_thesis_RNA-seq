import re

if __name__ == "__main__":

	Annotation_GO = "/home/binliu/R/Project_5_GOterms_Visualization/website_version/Annotation/mgi.gaf"
	All_GO = "/home/binliu/R/Project_5_GOterms_Visualization/website_version/Annotation/BP_related_GOterms.txt"
	output_list = "/home/binliu/R/Project_5_GOterms_Visualization/website_version/Annotation/All_GO_list.txt"
	output = open(output_list,"w")

	GO_node = []
	for line in open(All_GO, "r"):
		nodes = str(line).split()
		GO_child = nodes[0]
		GO_parent = nodes[2]
		if GO_child in GO_node: 
			continue
		else:
			GO_node.append(GO_child)
		if GO_parent in GO_node: 
			continue
		else:
			GO_node.append(GO_parent)
	#generate a GO_node array only include biological_process

	for line in open(Annotation_GO,"r"):
		node = str(line).split()
		Anno_GO = node[3]
		Anno_Gene = node[2]
		if Anno_GO in GO_node:
			outLine = Anno_GO + "\t" + Anno_Gene + "\n"
			output.write(outLine)
		else:
			continue			
	output.close()
	
