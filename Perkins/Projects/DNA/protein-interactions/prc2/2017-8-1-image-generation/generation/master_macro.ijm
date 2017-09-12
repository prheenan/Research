var global_iter_num = 0;

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

function get_line_id(line_array_x,line_array_y)
{
	// Args:
	//	line_array_<x,y>: see draw_array
	// Returns: 
	//	an id with the line given by x and y, assuming unique endpoints
	x0_y0=  x_y_str(line_array_x[0],line_array_y[0]);
	n = lengthOf(line_array_x);
	xn_yn = x_y_str(line_array_x[n-1],line_array_y[n-1]);
	return  "_" + x0_y0 + "_to_" + xn_yn; 
}

function save_and_draw_current_selection(suffix)
{
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
	output_dir = getDirectory("Image")
	getSelectionCoordinates(line_array_x,line_array_y);
	str_id =  get_line_id(line_array_x,line_array_y);
	output_path = output_dir+ getTitle() + str_id + suffix;
	print(output_path)
	saveAs("XY Coordinates", "" + output_path);
	draw_array(line_array_x,line_array_y);
}


macro "Save line segment as a dna-protein complex [C]"
{
	save_and_draw_current_selection("_Protein_DNA_complex_" + to_int_str(global_iter_num)); 
	global_iter_num += 1;
}
macro "Save line segment as a dna molecule [D]"
{
	save_and_draw_current_selection("_only___DNA_complex"+ to_int_str(global_iter_num));
	global_iter_num += 1;
}

