Please save or rename the GO term result generated from Panther with the name "Panther_result.txt"
The Panther_result document should be in the same folder with the Automation_visualization.sh document
Please enter "overrepresentation" or "enrichment" based on the analysis being run in Panther
The essential tables and documents for the generation of html webpage and the cytoscape division are automatically saved in the "result" folder
Before generating new results, it is recommended to delete the older ones to avoid confusion. But please do NOT delete the result folder! Otherwise it will cause unexpected errors. 

The documents generated in the "result" folder include:
  
additional_table.js - The result table in the JSON format. Would be used to generate a searchable table with essential information. 
 
best_pathway_nodes.sif - The .sif document which can be used in Cytoscape plugin to generate the network, the format should be adjusted manuelly. Recommended is an adjustment based on the "Default style". 

GOnrs_and_names_involved_in_best_pathway.txt - All the GO IDs and names included in the best pathway (including the ones being added to retain the connection between remote parent-child relationships). Could be integrated into the network by clicking Import - Table - File, please do not forget to click on the "Advanced Options.." and choose not to "Use first line as column names".

Panther_result_as_table.txt - The information table of the Panther result. Could also be integrated into the network by clicking Import - Table - File. The first line could be automatically used as column names. But please click on the "Advanced Options.." and choose the right delimiter "; semicolon" to make sure the table can be separated rightly.
 

