function tsvToTable(data){
	var lines = data.split("\n");
	var tbl = document.createElement('table');
	tbl.className = "result_table";
	tbl.border="1"
	var thead = document.createElement('thead');
	var tr = document.createElement('tr');
	tr.className = "evenrowcolor";
	var cells = lines[0].split("\t");
	var cdr3column = [];
	for(var a = 0;a < cells.length;++a){
		if(cells[a] == "CDR3 Sequence" || cells[a] == "CDR3_Sense_Sequence" || cells[a] == "Clone Sequence"){
			cdr3column.push(a);
		}
		var td = document.createElement('td');
		td.appendChild(document.createTextNode(cells[a]));
		tr.appendChild(td);
	}
	thead.appendChild(tr);
	tbl.appendChild(thead);
	var tbdy = document.createElement('tbody');
	
	for(var a = 1;a < lines.length;++a){
		tr = document.createElement('tr');
		var cells = lines[a].split("\t");
		if(cells.length == 1){
			continue;
		}
		for(var b = 0;b < cells.length;++b){
			td = document.createElement('td');
			td.appendChild(document.createTextNode(cells[b]));
			if(cdr3column.indexOf(b) != -1){
				td.className = td.className + " cdr3sequence"
			}
			tr.appendChild(td)
		}
		
		if(a % 2 == 0){
			tr.className = "evenrowcolor";
		} else {
			tr.className = "oddrowcolor";
		}
		
		tbdy.appendChild(tr);
	}
	tbl.appendChild(tbdy);
	return tbl;
}

function loadfile(file, patient, type){
	patient = patient.replace(".", "\\.");
	$('#hidden_div').load(file, function(){
		$('#result_div_' + patient + '_' + type).html(tsvToTable($('#hidden_div').html()));
		$('#result_div_' + patient + '_' + type + ' tr').hover(function() {
			$(this).addClass('hover');
		}, function() {
			$(this).removeClass('hover');
		});
		$('#result_div_' + patient + '_' + type + ' table').addClass('result_table');
		//$('#result_div_' + patient + ' tr:odd').addClass("oddrowcolor");
		//$('#result_div_' + patient + ' tr:even').addClass("evenrowcolor");
		$('#result_div_' + patient + '_' + type + ' table').before( "<a href='" + file + "'>Download " + file.replace(".txt", "") + "</a>" );
	});
}

var currentTD = new Array();

$( document ).ready(function() {
	$('.summary_table tr').hover(function() {
		$(this).addClass('hover');
	}, function() {
		$(this).removeClass('hover');
	});
	
	$('.summary_table tr:odd').addClass("oddrowcolor");
	$('.summary_table tr:even').addClass("evenrowcolor");
	
	$('.summary_table td[data-patient]').click(function() {
		var tmp = $(this);
		if(currentTD[tmp.attr("data-patient")] != null){
			currentTD[tmp.attr("data-patient")].removeClass("clicked_summary");
		}
		currentTD[tmp.attr("data-patient")] = tmp;
		currentTD[tmp.attr("data-patient")].addClass("clicked_summary");
	});
});
