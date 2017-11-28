
var global_iter_num = 0;
var global_protein_x = newArray;
var global_protein_y = newArray;
var global_protein_i = newArray;

function append(arr, value) 
{	
	// appends <value> to newArray <arr>. 
	// compatible with new and old-style imageJ
	arr2 = newArray(arr.length+1);
     	for (i=0; i<arr.length; i++)
	{
        		arr2[i] = arr[i];
	}
	arr2[arr.length] = value;
 	return arr2;
}

function clear_protein_arrays()
{
	// Removes all elements from the protein arrays
	global_protein_x = newArray();
	global_protein_y = newArray();
	global_protein_i= newArray();
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
		setForegroundColor(255,255,255);
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

function write_coords_to_file(output_path,line_array_x,line_array_y,header)
{
	// writes the current coodinates to the <output_path>
	// read the file for output
	file_var = File.open(output_path);
	// Get the total string to use
	newline = "\n";
	delim = " , " ;
	output_string = header + newline;
	n = lengthOf(line_array_x);
	for (i=0 ; i< n; i++)
	{
		tmp_x = "" + to_int_str(line_array_x[i]);
		tmp_y = "" + to_int_str(line_array_y[i]);
		output_string += "" + tmp_x + delim + tmp_y + newline;
	}
	output_string += "# (C) 2017 Patrick Heenan";
	// save the string to the output file
	print(output_string);
	File.close(file_var);
}

function save_and_draw_current_selection(suffix)
{
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
	output_dir = getDirectory("Image")
	getSelectionCoordinates(line_array_x,line_array_y);
	// already have a number for the ID 
	str_id =  ""
	output_path = output_dir+ getTitle() + str_id + suffix + ".txt";
	// Get the header of the protein locations (if any)
	header = get_protein_location_header();
	// write the formatted file...
	write_coords_to_file(output_path,line_array_x,line_array_y,header);
	// Remove any tagged locations
	clear_protein_arrays();
	draw_array(line_array_x,line_array_y);	
}

function add_current_point_as_protein()
{
	// Assuming we have a current selection, appends
	// the current x,y, and index associated with the selection
	// to the global arrays
	getSelectionCoordinates(line_array_x,line_array_y);
	current_length = lengthOf(line_array_x);
	// get the values we want
	current_i = current_length -1;
	current_x = line_array_x[current_i];
	current_y = line_array_y[current_i];
	// add them to the arrays
	global_protein_x = append(global_protein_x,current_x);
	global_protein_y = append(global_protein_y,current_y);
	global_protein_i = append(global_protein_i,current_i);
}

function get_protein_location_header()
{
	// Returns: The header desired for the protein locations
	n = global_protein_i.length;
	header = "# Protein locations <idx,x,y>: ";
	if (n == 0)
	{
		header += "No Proteins found.";
	}
	delim = ",";
	for (i=0; i<n; i++)
	{
		tmp_i =  to_int_str(global_protein_i[i]);
		tmp_x =  to_int_str(global_protein_x[i]);
		tmp_y =  to_int_str(global_protein_y[i]);
		header += "<" + tmp_i + delim  + tmp_x + delim  + tmp_y +">" + delim;
	}
	return header;
}

macro "Mark protein location [P]"
{
	add_current_point_as_protein();
	header = get_protein_location_header();
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


macro "Save line segment as a uninterpretable [U]"
{
	save_and_draw_current_selection("_Unknown_"+ to_int_str(global_iter_num));
	global_iter_num += 1;
}

