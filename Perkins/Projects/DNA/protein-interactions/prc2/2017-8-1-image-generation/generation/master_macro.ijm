
var global_iter_num = 0;
var global_protein_x = newArray(0);
var global_protein_y = newArray(0);
var global_protein_idx = newArray(0);

add_current_point_to_protein_arrays();

function add_current_point_to_protein_arrays()
{
	// Assuming we have a current selection, appends
	// the current x,y, and index associated with the selection
	// to the global arrays
	getSelectionCoordinates(line_array_x,line_array_y);
	current_length = lengthOf(line_array_x);
	// get the values we want
	current_idx = current_length -1;
	current_x = line_array_x[current_idx];
	current_y = line_array_y[current_idx];
	// add them to the arrays
	global_protein_x = Array.concat(global_protein_x,current_x);
	global_protein_y = Array.concat(global_protein_y,current_y);
	global_protein_i = Array.concat(global_protein_idx,current_idx);
}

function clear_protein_arrays()
{
	// Removes all elements from the protein arrays
	global_protein_x = Array.trim(global_protein_x,0);
	global_protein_y = Array.trim(global_protein_y,0);
	global_protein_idx = Array.trim(global_protein_idx,0);
}

function draw_array(line_array_x,line_array_y)
{
	// permanently draws a line on the image
	// Args:
	//	line_array_<x,y>: arrays, size N, of the <x,y> coords>
	// Returns:
	//	Nothing, but draws the line ontop of the image 
	for (i=0; i<lengthOf(line_array_x)-1; i++)
	{
		makeLine(line_array_x[i],line_array_y[i],line_array_x[i+1],line_array_y[i+1]);
		setForegroundColor(255,0,0);
		run("Draw");
	}
}

function to_int_str(num)
{
	// Returns: <num>, formatted as an integer string 
	return d2s(num,0);
}

function x_y_str(x,y)
{
	// Returns: the literal "_<x>_<y>", as a string (first thing must be a string, or it is interpretted as a number..
	ret_str = to_int_str(x);
	ret_str +=  "_" + to_int_str(y);
	return ret_str;
}


function save_and_draw_current_selection(suffix)
{
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
	output_dir = getDirectory("Image")
	getSelectionCoordinates(line_array_x,line_array_y);
	// already have a number for the ID 
	str_id =  ""
	output_path = output_dir+ getTitle() + str_id + suffix;
	print(output_path)
	saveAs("XY Coordinates", "" + output_path);
	draw_array(line_array_x,line_array_y);	
	// Remove any tagged locations
	clear_protein_arrays()
}


macro "Save line segment as a dna-protein complex [C]"
{
	save_and_draw_current_selection("_Protein_DNA" + to_int_str(global_iter_num)); 
	global_iter_num += 1;
}
macro "Save line segment as a dna molecule [D]"
{
	save_and_draw_current_selection("_DNA"+ to_int_str(global_iter_num));
	global_iter_num += 1;
}

macro "Save line segment as a multimer [M]"
{
	save_and_draw_current_selection("_Multiple_"+ to_int_str(global_iter_num));
	global_iter_num += 1;
}

macro "Mark protein location [P]"
{
	add_current_point_to_protein_arrays();
	n = lengthOf(global_protein_idx)
	print(n)
	for (i=0; i<n; i++)
	{
		print(global_protein_idx[i])
	}	
}

macro "Save line segment as a uninterpretable [U]"
{
	save_and_draw_current_selection("_Unknown_"+ to_int_str(global_iter_num));
	global_iter_num += 1;
}

