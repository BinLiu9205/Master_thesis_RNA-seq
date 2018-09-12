import re
import os


if __name__ == "__main__":
	cwd = os.getcwd()
	os.chdir(cwd)
	GOterm_library_single = {}
	GOterm_library = {}
	GO_id = "initialize"
	GO_synonym = []
	GO_relation = []
	GO_xref = []
	GO_other_relation = []
	GO_subset = []

	inputOBO = "go.obo"
	inputnrs = os.path.join(cwd,"result/all_nodes_in the pathway.txt")
	outputnames = os.path.join(cwd,"result/GOnrs_and_names_involved_in_best_pathway.txt")

	outname = open(outputnames, "w")
	
	for line in open(inputOBO, "r"):
		if line == "[Typedef]\n":
			break 
		if line == "[Term]\n":
			if GOterm_library_single != {}:
				GOterm_library[GOterm_library_single["id"]]=GOterm_library_single
				GOterm_library_single = {}
				GO_synonym = []
				GO_relation = []
				GO_xref =[]
				GO_other_relation = []
				GO_subset = []
		else:
			if str(line)[0:3] == "id:":   # would be the lowest level of dictionaries, save different aspects of single GO
				GO_id = str(line)[4:]	
				GOterm_library_single["id"] = GO_id
			if str(line)[0:4] == "def:":
				GO_def = str(line)[5:]
				GOterm_library_single["def"] = GO_def
			if str(line)[0:5] == "name:" :
				GO_name = str(line)[6:]
				GOterm_library_single["name"] = GO_name
			if str(line)[0:10] == "namespace:" : 
				GO_cate = str(line)[11:]
				GOterm_library_single["category"] = GO_cate
			if str(line)[0:8] == "synonym:" :
				GO_synonym_s = str(line)[9:]
				GO_synonym.append(GO_synonym_s)
				GOterm_library_single["synonym"] = GO_synonym
			if re.search("^is_a:",line) != None: 
				GO_relation_s = str(line)
				GO_relation.append(GO_relation_s)
				GOterm_library_single["relationship"] = GO_relation
			if re.search("^intersection_of:",line) != None or re.search("^relationship:",line) != None: 
				GO_other_relation_s = str(line)
				GO_other_relation.append(GO_other_relation_s)
				GOterm_library_single["relationship_additional"] = GO_other_relation
			if str(line)[0:5] == "xref:" :
				GO_xref_s = str(line)
				GO_xref.append(GO_xref_s)
				GOterm_library_single["xref"] = GO_xref
			if str(line)[0:7] == "subset:" :
				GO_subset_s = str(line)[8:]
				GO_subset.append(GO_subset_s)
				GOterm_library_single["subset_of"] = GO_subset

	# add the last one:			
	GOterm_library[GOterm_library_single["id"]]=GOterm_library_single

	for line in open(inputnrs):
		outLine = str(line)[0:10] + "\t" + GOterm_library[str(line)]["name"] 
		outname.write(outLine)
	
	#print(GOterm_library)

	
	

	


	
