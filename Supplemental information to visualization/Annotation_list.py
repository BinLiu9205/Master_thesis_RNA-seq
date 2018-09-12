import re

if __name__ == "__main__":
	input_Anno  = "/home/binliu/R/Project_5_GOterms_Visualization/website_version/Annotation/mgi.gaf"
	input_BP = "/home/binliu/R/Project_5_GOterms_Visualization/BP_related_GOterms.txt"
	output_Anno = "/home/binliu/R/Project_5_GOterms_Visualization/website_version/Annotation/full_annotation.txt"
	Anno_in = open(input_Anno, "r") 
	BP_in = open(input_BP, "r") 
	Anno_out = open(output_Anno, "w")
	Annotation = {} 
#save the GO_gene relationship into memory
	for line in Anno_in:
		nodes = str(line).split()
		GO_id = nodes[3]
		Gene_sym = nodes[2]
		if GO_id not in Annotation:
			Annotation[GO_id] = []
			Annotation[GO_id].append(Gene_sym)
		else:
			Annotation[GO_id].append(Gene_sym)

#find direct children parent relationship to each other
#parent to children, key-value pair

	relation = {}
	for line in BP_in: 
		nodes = str(line).split()
		child = nodes[0]
		parent = nodes[2]
		if child not in relation:
			relation[child] = []
			relation[child].append(parent)
		else:
			relation[child].append(parent)

#expand direct relationship to remote relationship

	for c_p in relation: 
		dir_list = relation[c_p]
		for i in range(0,len(dir_list)):
			r_node = dir_list[i]
			if not (r_node) == 'GO:0008150':
				new_relationship = relation[r_node]
				for j in range(0,len(new_relationship)):
					new_r_node = new_relationship[j]
					if new_r_node not in relation[c_p]:
						relation[c_p].append(new_r_node)
					else: continue
			else:
				continue

#child to all its parents, including remote parents, add annotation of child back to his parents
	for c_c in relation:
		par_list = relation[c_c]  # find out all the parents of single child 
		if c_c in Annotation:
			save_anno = Annotation[c_c] # save the annotation of single child into save_anno
			for i in range(0,len(par_list)):
				r_node = par_list[i] # parent being splitt, and run seperately
				if r_node in Annotation:
					for j in range(0,len(save_anno)): 
						node_anno = save_anno[j] #every annotation of child being saved
						if node_anno not in Annotation[r_node]:
							Annotation[r_node].append(node_anno)
						else: continue
				else: continue
		else:
			continue
	


	
			

    

