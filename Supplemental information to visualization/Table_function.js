var globalCounter = 0;
var mytable = document.getElementById('mytable');
for (var i = 0; i < obj.length; i++) {
	var tr = "<tr";

	tr += " id=" +  obj[i].GO_ID + ">" +
	"<td>" + obj[i].GO_ID + "</td>" + "<td>" +obj[i].GO_Name  + "</td>" + "<td>" + obj[i].Size.toString() + "</td>"+ 
	"<td>" + obj[i].DE_Genes.toString() + "</td>" + "<td>" + obj[i].Expected.toString() + "</td>" + 
	"<td>" + obj[i].Direction + "</td>" + "<td>" + obj[i].Enrichment.toString() + "</td>" + 
	"<td>" + obj[i].p_value.toString() + "</td></tr>";

	/* We add the table row to the table body */
	mytable.innerHTML += tr;
} ;
   
var previousValue = " ";
$('tr').click( function (e) {
 cy.nodes("[shared_name='"+previousValue+"']").unselect();
 previousValue = $(this).find('td:first').text();
 //alert(previousValue);
 var n = cy.nodes("[shared_name='"+previousValue+"']");

 //n.style("background-color" , "rgb(255,255,0)") // this just makes it yellow

 $("body").scrollTop(0);
 n.select();
 //alert(cy.nodes(n));

});   	

