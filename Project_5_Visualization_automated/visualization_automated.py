import re
import os 

if __name__ == "__main__":
	cwd = os.getcwd()
	os.chdir(cwd)
	inputOBO = os.path.join(cwd,"go.obo")
	outputBP_term = os.path.join(cwd,"result/BP_related_GOterms.txt")
	out_terms = open(outputBP_term,"w")  #make them able to be written and modified
	
	newTerm = False
	GOvalue = "initilize"
	isAfound = False
	isABP = False

	for line in open(inputOBO, "r"):
		if newTerm:
			GOvalue = line[4:-1]
			newTerm = False
			isAfound = False 
			isABP = False

		if line == "[Term]\n":
			newTerm = True
		if line == "[Typedef]\n":
			break   #read all the nodes, don't need type defs
		if "biological_process" in line:
			isABP = True
		if line == "\n":    #end of entry
			isAfound = False
		elif re.search("^is_a:",line) != None:
			isAfound = True
        
		if isABP:
			if isAfound:
				splits = line.split()
				outLine="\n"
				if splits[0] == "is_a:":  
					outLine = GOvalue + "\t" + splits[0] + "\t" + splits[1] + "\n"
					out_terms.write( outLine )
	out_terms.close()
             
##above was about producing a GOterms list which only include Biological Process, only is_a relationship would be selected for further analysis###

    
	inputBP_term = os.path.join(cwd,"result/BP_related_GOterms.txt")
	output_relation = os.path.join(cwd,"result/parent_child.txt")
	out_relation = open(output_relation, "w")
   
	newNode = True
	old_child = "initialize"   
	s_parent = "initialize"
	child = "initialize"
	parent = []
	GO_term_graph = {}
    

	def merge_dict (x, y):
		x = x.copy()
		x.update(y)
		return x

	infile = open(inputBP_term, "r")
	for line in infile:
		nodes = str(line).split()# split into a list [child, , parent]   
		child = nodes[0]
		s_parent = nodes[2]
		newNode = old_child != child
		if newNode:
			GO_term_graph[old_child] = parent # save the previous child to the dictionary
			outNode = old_child + "\t" + str(parent) + "\n"
			out_relation.write(outNode)
			parent = [s_parent]
			old_child = child  
		else:
			parent.append(s_parent)
	GO_term_graph[child] = parent
	outNode = old_child + "\t" + str(parent) + "\n"
	out_relation.write(outNode)
    
	out_relation.close()      
    #print("Counter:",counter)
        #if not old_child == child:
           #newNode = True
        #else:
           #newNode = False
    
	#print (GO_term_graph)
    #print ("---------------This is the end of GO_term_graph construction part)----------------")

##### Codelines above were to save GOterms in the memory, the parent of every child would be save in a list corresponding ####

	def find_all_paths(graph, start, end, path=[]):
		path = path + [start]
		if start == end:
			return [path]
		if not start in graph:
			return []
		paths = []
		for node in graph[start]:
			if node not in path:
				newpaths = find_all_paths(graph, node, end, path)
				for newpath in newpaths:
					paths.append(newpath)
		return paths
   
    #print (find_all_paths(GO_term_graph, "GO:2001307", "GO:0044249")) #be used to test whether the definition was right

    #test = "/home/binliu/R/Project_5_GOterms_Visualization/test.txt"
    #test_file = open(test, "w")
    #test_file.write = find_all_paths(GO_term_graph, "GO:2001317", "GO:2001316")

   
### def was a code based on https://www.python.org/doc/essays/graphs/ , it is possible to find all the list connecting two nodes ####
### could be used for finding common ancester ?? #########

        

##GO:0008150 would be the GO_id of biological_process, theoritically, all the GOterms related in this area would be somehow connected with it

	GO_sig = os.path.join(cwd,"result/GO_number_for_visualization.txt") # the significant GO_terms list
	#print(os.path.join(cwd,"GO_number_for_visualization.txt")) 
	GO_path = os.path.join(cwd,"result/GO_pathway.txt") 
	best_path = os.path.join(cwd,"result/GO_best_pathway_for_cytoscape.sif") 
	Nr_for_name = os.path.join(cwd,"result/GO_numbers_involved_in_best_pathway.txt") 

	out_path = open(GO_path, "w")
	out_best = open(best_path, "w") 
	for_name = open(Nr_for_name, "w")
	char_GO_sig = []
	pathway_related = []
	pathway_sub = []
	all_path_sub = {}
	counter = 0
	max_counter = -1
	sig_path = 0
	last_node = "initialize"
	outLine = "/n"
	outLine_2 = "/n"

	for line in open(GO_sig, "r"):
		single_line = str(line).split()
		single_GO = single_line
		char_GO_sig.append(str(str(single_GO)[2:12])) # transfer the GOterms into characters and save them into a list 
        #print(len(char_GO_sig))

	for GO in range(0,(len(char_GO_sig))):
		pathway_related = find_all_paths(GO_term_graph, char_GO_sig[GO] , 'GO:0008150')
		if len(pathway_related) ==1:                       # when there is only one pathway, loop will out of range, write down this path
			pathway_sub = []
			pathway_sub.append(pathway_related[0])
			out_path.write(char_GO_sig[GO] + "\t"+ str(pathway_sub) + "\n") 
			all_path_sub[char_GO_sig[GO]] = pathway_sub[0]
	
		else:
			sig_path = 0
			for i in range(0,len(pathway_related)):    # how many levels are between them when reaching biological_process
				this_path = pathway_related[i]
				if len(this_path) < 3: print("Warning: pathway len was less than 3:" , this_path)
                  #count the number of significant GO terms in this_pathway
				for j in range(1,(len(this_path)-1)):
					if this_path[j] in char_GO_sig:
						counter += 1
				if counter > max_counter:
					max_counter= counter
					sig_path = i
				if sig_path == 0: print("WARNING: sig_path might be default ",GO,char_GO_sig[GO],pathway_related)
				best_path = pathway_related[sig_path]
				out_path.write(char_GO_sig[GO] + "\t"+ str(best_path) + "\n")  
              # finally, add this best pathway for this go term to a dictionary:
				all_path_sub[char_GO_sig[GO]] = best_path                     
	out_path.close()
     
     
    #print(all_path_sub)
    #print(len(all_path_sub))

	for GO in range(0,(len(char_GO_sig))): # for every go term, get the best pathway from the dictionary
		path_value = all_path_sub[char_GO_sig[GO]]
        #print(char_GO_sig[GO])
        #if len(path_value)==2:
                 #outLine = path_value[0] + "\t" + path_value[1] + "\n"
                 #out_best.write(outLine)
        #else
		if len(path_value) < 2 : print ("Warning: path_value was less than 2")
		for nr in range(0,len(path_value)): # for every term in this pathway
			outLine_2 = path_value[nr] + "\n"
			for_name.write(outLine_2)
			if nr == 0:
				last_node = path_value[0]
			else:
				outLine = last_node + "\t" + "is_a:" + "\t" + path_value[nr] + "\n" 
				out_best.write(outLine)
				last_node = path_value[nr]      
	out_best.close()
	for_name.close()
            


