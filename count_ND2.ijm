// Cell counting program 
// - - - - - - - - - - 

/*
   COPYRIGHT DISCALIMER:
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.



   Author: Christian Evenhuis, christian.evenhuis@gmail.com
   Date:   August 2017
*/

// Update the macros, remember to:
//    Help > Refresh Menus
//str=exec("/Users/evenhuis/Dropbox/UTS/imageJ_macros/assemble_macros.sh");

//- - - Global parameter with their defaults- - - - - - 
var file_ext='nd2';
var px2um=NaN;
var width_um, height_um;		// size of the image in um

var process_chan=1;		// channel to process thresholds on

var nseq=0;				// Number of series in ND2 file
var nchan=1;
var save_overlay=true;	// Save a summary slide
var nopen=1;			// series number to open
var thresh_upper=2.5;
var watershed=false;
var size_min="0";
var size_max="Infinity";
var circ_min=0.0;
var circ_max=1.0;


var use_fit=false;
var Rhalf_max="Infinity";
var R2_min=0.0;
var depth;             // depth of sample

var repeat=true; // loop variable for parameter setting
var Im_x0=10;       // location of the image
var Im_y0=100; 
var Im_w0=800; 
var Im_h0=800;

var dil_fac=1.0;	// dilution applied

var nres;	// number of cells found

var pi=3.141592654;

var slash=File.separator;  // Stupid windows!
// - - -  - end Global definition - - -  - -  - -  - - - -- - - -- 

//process_chan=3;
//id=getImageID()
//threshold_image(id,-4.);

//calibrate_cell_count_f();
//count_ND2_f();
process_directory();
exit;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
macro "calibrate_cell_count [C]" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	setBatchMode(false);
	calibrate_cell_count_f();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
macro "perform_cell_count [P]" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	setBatchMode(false);
	count_ND2_f();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
macro "process directory [D]" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	setBatchMode(false);
	process_directory();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function process_directory() {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path 	=getDirectory("Select a directory");	// Name of file 
print(path);
succ = look_for_previous_parameters(path);

// get files in directory
print("file ext",file_ext);
files=list_files_in_dir(path,file_ext);
for( i=0 ; i<files.length ; i++ ){
		
	file = files[i];
	
	print(file);
	filepath=path+file;
	count_file(filepath);
}

return;
} // end process dir


// - - - - - - - - - - - - - - - - - - - - - - - -
function fit_edge( id, n0, n1 ){
// - - - - - - - - - - - - - - - - - - - - - - - -

selectImage(id);
Stack.setPosition(process_chan,1,1);

nmax=pow(2,bitDepth());
if(nmax==256){
	getHistogram(vals,counts,256);
}else{
	getHistogram(vals,counts,256,0,nmax);
}
tmp=get_hist_mode_lhm_hhm(counts);


//if( tmp[2]==256 ){
	//setBatchMode(false);
	//exit;
//}
mode=vals[minOf(maxOf(tmp[0],0),255)]; 
fwhm =vals[minOf(maxOf(tmp[2],0),255)]-vals[minOf(maxOf(tmp[1],0),255)];

dist=newArray(12);
grey=newArray(12);
for( i=n0; i<n1; i++ ){
	x=getResult("X",i)/px2um;
	y=getResult("Y",i)/px2um;
	major =  getResult("Major",i)/2./px2um;
	minor =  getResult("Minor",i)/2./px2um;
	angle = -getResult("Angle",i)*pi/180.;

	rmin=0.9*major;
	rmax=2.5*major;
	dr=(rmax-rmin)/11.;


	run("Select None");
	//setBatchMode(false);
	for( k=0 ; k<12; k++ ){
		r1=maxOf(rmin + (k-1)*dr,0);
		r2=  r1 + dr;
		
		// make an annular selection
		makeEllipse(x+r2*cos(angle),y+r2*sin(angle),
		            x-r2*cos(angle),y-r2*sin(angle),
		            minor/major);
	    setKeyDown("alt");
		makeEllipse(x+r1*cos(angle),y+r1*sin(angle),
		            x-r1*cos(angle),y-r1*sin(angle),
		            minor/major);	
		dist[k]=r1-rmin;

		grey[k]=getMedian();
		//print(dist[k]+" "+grey[k]);
		//waitForUser("d");
	}
	run("Select None");

	// work out the half way point
	Imin=grey[0];
	Ihalf=0.5*(mode+Imin);
	khalf=11;
	for(k=1;k<12; k++ ){
		if( grey[k]>Ihalf ){
			khalf=k;
			k=12;
		}
	}
	//print("Rhalf=",dist[khalf]*px2um,dist[0]*px2um);
	setResult("cell_Dark",i,Imin);
	setResult("R_half",i,(dist[khalf]-dist[0])*px2um);
	//updateResults();
	//waitForUser("d");
}


updateResults();
return;
	
}

// - - - - - - - - - - - - - - - - - - - - - - - -
function getMedian(){
// - - - - - - - - - - - - - - - - - - - - - - - -	
	getStatistics( area, mean, min, max, std, hist );
	// look for the median intensity (this is the area/2 sort pixel value
	med=area/2;
	cdf=0;
	imed=0;
	for( k=0; k<255; k++ ){
		if( cdf<= med && med<=cdf+hist[k] ){
			imed=k;
			k=256;
		}else{
			cdf+=hist[k];
		}
	}
	return(imed);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function count_ND2_f(){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Get the user to choose a file
// Get the ND2 file to open
filepath=File.openDialog("Select a File");	// Name of file 
count_file(filepath);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function count_file( filepath){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path    =File.getParent(filepath);	     	// directory it is in
file    =substring(filepath,lengthOf(path)+1);
ne      =lastIndexOf(file,".");
name    =substring( file,0,ne);


cfile=path+slash+name+".cell_counts.txt";
if( File.exists(cfile)){
	print(name+".cell_counts.txt : already exits!");
	return;
}

initial_setup(filepath);
look_for_previous_parameters(path);

dil_fac = get_dilution_factor_from_filename(file);


cell2conc=1E-3     * 1E-3  /((width_um-2*size_max)*1E-6*(height_um-2*size_max)*1E-6*depth*1E-6);
print("celAls to cell/mL ",cell2conc);



// processing loop
setBatchMode(true);
for(ns=1;ns<=nseq;ns++){
	open_series( filepath, ns );
	w_id=getImageID();
	setLocation(Im_x0, Im_y0, Im_w0, Im_h0);


	t_id = threshold_image       ( w_id, thresh_upper );
	//t_id= threshold_image_local_loc( w_id, -thresh_upper, 1 );

	binary_process_image( t_id );

	n0=nResults();
	measure_cells_in_image( t_id, ns );
	n1=nResults();
	
	selectImage(w_id);
	fit_edge     (w_id,n0,n1);
	calculate_RGB(w_id,n0,n1);
	
	// add filters
	for(i=n0;i<n1;i++){
		setResult("Exclude",i, apply_filters(i));
	}
	updateResults();
	roiManager("reset");
	
	selectImage(t_id ); 	
	close();
	selectImage(w_id);

	add_count_summary(ns,n0,n1);
	display_counted_cells(n0,n1);
	
	// Update the progress bar
	js=2;
	ps = ns%js;
	if( ps==0 ){
		showProgress((i-js*ps)/nseq);
		run("Collect Garbage");
	}
}
selectImage(w_id);
close();
setBatchMode(false);

saveAs("Results",path+slash+name+".cell_sizes.txt");
IJ.renameResults("Results","Sizes");

IJ.renameResults("Counts", "Results");
saveAs("Results",path+slash+name+".cell_counts.txt");
IJ.renameResults("Results","Counts");


}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function calibrate_cell_count_f(){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Get the user to choose a file
// Get the ND2 file to open
filepath=File.openDialog("Select a File");	// Name of file 
path    =File.getParent(filepath);	     	// directory it is in
file    =substring(filepath,lengthOf(path)+1);


initial_setup(filepath);
look_for_previous_parameters(path);
get_user_parameters(" ");

// Open a series
//setBatchMode(true);
open_series( filepath, nopen );
w_id=getImageID();
setLocation(Im_x0, Im_y0, Im_w0, Im_h0);


// processing loop
repeat=true;
nopen_old=nopen;
loop=0;
while( repeat ){
	nopen=minOf(nopen,nseq);
	nopen=maxOf(nopen,1);
	if( ! (nopen_old == nopen) ){
		selectImage(w_id); close();
		setBatchMode(false);
		open_series( filepath, nopen );
		w_id=getImageID();
		//setBatchMode(false);
		setLocation(Im_x0, Im_y0, Im_w0, Im_h0);
	}

	selectImage(w_id);
	Overlay.clear();
	
	setBatchMode(true);
	selectImage(w_id);
	run("Duplicate...", "duplicate");
	d_id = getImageID();
	

	t_id = threshold_image( d_id, thresh_upper );

	binary_process_image( t_id );


	run("Clear Results");
	n0=nResults();
	measure_cells_in_image( t_id, nopen );
	n1=nResults();
	

	fit_edge     (d_id,n0,n1);
	calculate_RGB(d_id,n0,n1);
	
	// add filters
	for(i=n0;i<n1;i++){
		setResult("Exclude",i, apply_filters(i));
		//setResult("Exclude",i,0);
	}
	updateResults();
	roiManager("reset");

	
	selectImage(t_id ); 	close();
	selectImage(d_id ); 	close();
	selectImage(w_id);
	setBatchMode(false);

	if(nchan>1){Stack.setChannel(process_chan)};
	add_count_summary(nopen,n0,n1);
	display_counted_cells(n0,n1);
	
	nopen_old=nopen; 	
	get_user_parameters(n1-n0);  // repeat is set to false in here
	nopen=minOf(nopen,nseq);
	nopen=maxOf(nopen,1);

	run("Collect Garbage");
	loop++;
}

file_ext=get_ext(filepath);
save_parameters(path);

}// end calibrate



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function add_count_summary(ns,n0,n1){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// add the number of cells to the count list
ntotal=n1-n0;
ncount=0; nsize=0;nsoft=0;nedge=0;
for(i=n0;i<n1;i++){
	if(getResult("Exclude",i)==0){
		ncount++;
		// calculate the volume
		maj=getResult("Major",i);
		min=getResult("Minor",i);
		setResult("Volume_um^3",i,4.*pi/3.*maj/2.*min/2.*(maj+min)/4.);
	};
	if(getResult("Exclude",i)==1){nsize --;};
	if(getResult("Exclude",i)==2){nsoft --;};
	if(getResult("Exclude",i)==3){nedge --;};
}
	
IJ.renameResults("Results","Sizes");
setLocation(Im_x0+Im_w0, Im_y0);

// Save the count results to a window


Ext.setSeries(ns-1);
Ext.getPlanePositionX(X0,0);
Ext.getPlanePositionY(Y0,0);

IJ.renameResults("Counts", "Results");
nr=nResults();
setLocation(Im_x0+2*Im_w0, Im_y0);
setResult("File_name",   nr,file);
setResult("Series",      nr,ns);
setResult("X0",          nr, X0);
setResult("Y0",          nr, Y0);
setResult("total",       nr,ntotal);
setResult("ex. size",    nr,nsize);
setResult("ex. soft",    nr,nsoft);
setResult("ex. edge",    nr,nedge);
setResult("counted",     nr,ncount);
print("depth "+depth);
setResult("conc cell/mL",nr,d2s(dil_fac*ncount*1E-3     * 1E-3  /((width_um-2*size_max)*1E-6*(height_um-2*size_max)*1E-6*depth*1E-6),-4));
								             // L->ml   m^3->L          W              um->M         H              um->  D mm-> M
setResult("dilution",    nr,dil_fac);
IJ.renameResults("Results", "Counts");
setLocation(Im_x0+Im_w0, Im_y0);
IJ.renameResults("Sizes","Results");
return;
}	// end add counts

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function initial_setup(filepath){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
run("Bio-Formats Macro Extensions");
Ext.setId(filepath);
Ext.getSeriesCount(seriesCount);


// Get the meta data
meta_str = get_meta_data( filepath );					// get the metadata string
nseq      = parseInt  (extract_key_val("uiSequenceCount",meta_str)); 	// work out the number of images
if(isNaN(nseq)){nseq=1;}
print("nseq",nseq);
nchan     = parseInt  (extract_key_val("SizeC",          meta_str)); 
px2um     = parseFloat(extract_key_val("dCalibration",   meta_str));  	// Conversion factor from um to pixel

if( isNaN(px2um)){
	
	path    =File.getParent(filepath);
	look_for_previous_parameters( path );
	print("after looking",px2um);
	if( isNaN(px2um) ){
		open_series( filepath, 1 );
		setTool("line");
		line_not_measured=true;
		while( line_not_measured){
			waitForUser("Draw a line along a grid edge");
			getSelectionCoordinates(xpoints, ypoints);
			if( xpoints.length==2){

				dx2 = (xpoints[1]-xpoints[0])*(xpoints[1]-xpoints[0]);
				dy2 = (ypoints[1]-ypoints[0])*(ypoints[1]-ypoints[0]);
				px2um = 250./sqrt( dx2+dy2);
			    line_not_measured=false;                  
			}                       
	}
}
width_um  = parseInt  (extract_key_val("SizeX",          meta_str))*px2um; 
height_um = parseInt  (extract_key_val("SizeY",          meta_str))*px2um;

size_max = minOf( width_um, height_um)*0.05;	// default max size (5% of image)
depth = 100.;		// haemo depth
//if( width_um > 800 ){ depth=1000.;}// guess the depth based on the size of the picture


// Close existing results windows that might ve lying around
test_and_close("Results");
test_and_close("Sizes");
test_and_close("Counts");
test_and_close("Summary");

// Set up the results tables
setResult("File_name",0," ");
IJ.deleteRows(0, 0); 
IJ.renameResults("Results","Counts");

run("Clear Results"); 
setResult("File_name",0," ");
IJ.deleteRows(0, 0); 
return;
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function calculate_RGB(id,n0,n1){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// add the RGB info
selectImage(id);

getDimensions(width, height, channels, slices, frames);



Stack.setPosition(1,1,1);
// Select all cells
run("Select None");
print(n1-n0);
if( n1-n0>0){
	for(i=0;i<n1-n0;i++){
		setKeyDown("shift");
		roiManager("select",i);
	}
	run("Make Inverse");
}else{
	run("Select All");
}
backgrounds=newArray(channels+1);
for(j=1;j<=channels;j++){
	Stack.setChannel(j);
	getStatistics(area, mean, min, max, std, histogram);
	nmax=pow(2,bitDepth());
	if(nmax==256){
		getHistogram(vals,counts,256);
	}else{
		getHistogram(vals,counts,256,0,nmax);
	}
	tmp=get_hist_mode_lhm_hhm(counts);
	mode=vals[minOf(maxOf(tmp[0],0),255)]; 
	backgrounds[j]=mode;
	print("back ground Chan"+j+" "+mean);
	//setBatchMode(false);
	//exit;
}

for(j=1;j<=channels;j++){
	Stack.setPosition(j,1,1);
	col="Chan"+j;	
	for(i=0;i<n1-n0;i++){
		roiManager("select",i);
		getStatistics(area,mean);
		setResult(col,         n0+i,round(mean));
		setResult(col+"_bgc",  n0+i,round(abs(mean-backgrounds[j])));
	}
}
run("Select None");
return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function apply_filters(i){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	// exclude on length
	major=parseFloat(getResult("Major",i));
	if(  major<size_min || size_max < major ){
		return(1);
	}

	// exclude on soft edge
	Rhalf=getResult("R_half",i);
	if(  Rhalf_max < Rhalf  ){ 
		return(2);
	}

	// Exclude those touching the edges
	// bounding box coordinates
	X0=getResult("BX",i);
	Y0=getResult("BY",i);
	X1=X0+getResult("Width",i);
	Y1=Y0+getResult("Height",i);

	// exclude those that are completely out
	if( Y1<size_max || height_um-size_max<Y0 || 
		X1<size_max ||  width_um-size_max<X0 ){ 
		return(3);	
	}
	
	// exlude the bottom completely
	if( Y0<height_um-size_max && height_um-size_max< Y1 ){
	 	count=3; 
		return(3);		
	}
	// exlude the left edge, minus the box at the top
	if( X0<size_max && size_max< X1 && size_min<Y0){
		return(3);		
	}
	
	return(0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function open_series( file_path, ns ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if( 0<=ns && ns<=nseq ){	
	bio_options="autoscale color_mode=Default concatenate_series  view=Hyperstack stack_order=XYCZT";
	run("Bio-Formats Importer", "open=["+filepath+"] "+bio_options+"  specify_range series_"+ns);
	setLocation(Im_x0, Im_y0, Im_w0, Im_h0);
	run("Set Scale...", "distance=1. known="+px2um+" unit=micron");
	getDimensions(width, height, channels, slices, frames);
	width_um = width*px2um;
	height_um = height*px2um;
}else{
	print("sequence out of range");
	exit;
}

return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function get_user_parameters(nfound){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
// create dialog box for options
nres=nResults();
Dialog.create("Cell counter and sizer");

Dialog.addNumber("Series number : ",nopen );
Dialog.addString("Number of cells found: ",nfound );
//Dialog.addRadioButtonGroup("Threshold", newArray("Auto","Manual"), 2, 1, "Manual");
Dialog.addNumber("Process channel ", process_chan ); 
Dialog.addNumber("Threshold", thresh_upper, 2, 5, "arb"); 

Dialog.addCheckbox("Watershed", watershed);
//var corner_correct=false;
//Dialog.addCheckbox("Correct corner", true)print(f,"Rhalf_max= "   +d2s(Rhalf_max,3));
//Dialog.addCheckbox("Count cells", true);
//Dialog.addCheckbox("Measure radius", true);
//Dialog.addCheckbox("Save a summary", save_overlay);


Dialog.addNumber("Mininum Length", size_min, 2, 10, "um");
size_min = maxOf(px2um*2,size_min); 
Dialog.addNumber("Maximum Length", size_max, 2, 10, "um"); 


Dialog.addNumber("Mininum Circ", circ_min, 2, 10, ""); 
Dialog.addNumber("Maximum Circ", circ_max, 2, 10, ""); 

//Dialog.addCheckbox("Use fit params (caution slow)", use_fit);
Dialog.addNumber("Softness of edge in um", Rhalf_max, 2, 10, ""); 
//Dialog.addNumber("Minimum R2", R2_min, 2, 10, ""); 

Dialog.addMessage("Haemocytomer   =  100um");
Dialog.addMessage("Segwick Rafter = 1000um");
Dialog.addNumber("Depth of sample (um)", depth, 0, 10, ""); 

Dialog.addCheckbox("I'm happy with the parameters", false)

Dialog.show();


// Process the options
//thresh_choice=Dialog.getRadioButton();
nopen       =Dialog.getNumber();
process_chan=Dialog.getNumber();
thresh_upper=Dialog.getNumber();

watershed      = Dialog.getCheckbox();
//corner_correct = Dialog.getCheckbox();
//count_cells    = Dialog.getCheckbox();
//measure_size   = Dialog.getCheckbox();
//save_overlay     =  Dialog.getCheckbox();

size_min = parseFloat(Dialog.getNumber());
size_max = parseFloat(Dialog.getNumber());

circ_min = parseFloat(Dialog.getNumber());
circ_max = parseFloat(Dialog.getNumber());

//use_fit		= Dialog.getCheckbox();
Rhalf_max	= parseFloat(Dialog.getNumber());
//R2_min		= Dialog.getNumber();

depth = parseFloat(Dialog.getNumber());

param_done   = Dialog.getCheckbox();
if( param_done ){ repeat=false; }	

size_min=maxOf(size_min,sqrt(16/pi)*px2um);

return(repeat);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function look_for_previous_parameters( path ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// Look for previous parameters in this file
param_file=path+slash+"cell_count.param.txt";
//print("param: "+param_file);

loaded=false;
if( File.exists(param_file) ){
	str=File.openAsString(param_file);
	lines=split(str,"\n");
	for( i=0; i<lines.length; i++ ){
		args=split(lines[i],"= ");

		//if( args[0]=="save_overlay"){ save_overlay=args[1]; }
		if( args[0]=="file_ext"    ){ file_ext=    args[1] ; }
		if( args[0]=="process_chan"){ process_chan=parseInt  (args[1]); }
		if( args[0]=="size_min"    ){ size_min=    parseFloat(args[1]); }
		if( args[0]=="size_max"    ){ size_max=    parseFloat(args[1]); }
		if( args[0]=="circ_min"    ){ circ_min=    parseFloat(args[1]); }
		if( args[0]=="circ_max"    ){ circ_max=    parseFloat(args[1]); }
		if( args[0]=="watershed"   ){ watershed=     args[1]; }
		//if( args[0]=="use_fit"     ){ use_fit=     args[1]; }
		if( args[0]=="Rhalf_max"   ){ Rhalf_max=   parseFloat(args[1]); }
		
		if( args[0]=="depth"   ){ depth=   parseFloat(args[1]); }
		if( args[0]=="px2um"   ){ px2um=   parseFloat(args[1]); print("set px2um",px2um); }
		//if( args[0]=="R2_min"      ){ R2_min=      args[1]; }
	}
	size_min=maxOf(size_min,sqrt(16/pi)*px2um);
	loaded=true;
}


return( loaded);
} // end look_for previous_parameters

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function save_parameters(path){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Save the parameters

param_file = path+slash+"cell_count.param.txt";
f=File.open(param_file);
//print(f,"save_overlay= " +save_overlay);
print(f,"file_ext= "+file_ext);
print(f,"process_chan= "+process_chan);
print(f,"thresh_upper= "+thresh_upper);
print(f,"size_min= "    +size_min);
print(f,"size_max= "    +size_max);
print(f,"circ_min= "    +d2s(circ_min,3));
print(f,"circ_max= "    +d2s(circ_max,3));
print(f,"watershed="    +watershed);
//print(f,"use_fit= "     +use_fit);
print(f,"Rhalf_max= "   +d2s(Rhalf_max,3));
print(f,"depth= "       +d2s(depth,0));
print(f,"px2um= "       +d2s(px2um,4));
//print(f,"R2_min= "      +d2s(R2_min,3));
File.close(f);

return;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function threshold_image( id, t_u ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Do the thresholding
selectImage(id);
Stack.setPosition(process_chan,1,1);

nmax=pow(2,bitDepth());
if(nmax==256){
	getHistogram(vals,counts,256);
}else{
	getHistogram(vals,counts,256,0,nmax);
}
tmp=get_hist_mode_lhm_hhm(counts);

mode=vals[minOf(maxOf(tmp[0],0),255)]; 
fwhm =vals[minOf(maxOf(tmp[2],0),255)]-vals[minOf(maxOf(tmp[1],0),255)];


run("Select All");
run("Copy");

run("Internal Clipboard");
t_id=getImageID();
run("Set Scale...", "distance=1. known="+px2um+" unit=micron");

if(t_u>0){
	setThreshold(0, mode-t_u*fwhm);
}else{
	setThreshold(mode-t_u*fwhm, nmax);
}
setOption("BlackBackground",false);
run("Convert to Mask");


return(t_id);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function binary_process_image( id ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// binary processing section

run("Fill Holes");
if(watershed){	
	run("Erode");
	run("Watershed");
	//run("Dilate");
}

return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function threshold_image_local_loc( id, t_l, ndiv ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Do the thresholding
selectImage(id);

bmode_on=is("BatchMode");

//Stack.setPosition(process_chan,1,1);

s_id=id;
getDimensions(w,h,ch,sl,fr);

selectImage(s_id); 

nx=ndiv;
ny=round(h/w*ndiv);
thresh_low=3.0;


Overlay.clear();

dx=w/nx;
dy=h/ny;


//if( ! bmode_on){ setBatchMode(true);};

selectImage(s_id);
run("Select All");
run("Copy");
run("Internal Clipboard");
st_id=getImageID();


newImage("Thesholded","16-bit",w,h,1);
tt_id=getImageID();

for( i=0 ; i<nx ; i++ ){
for( j=0 ; j<ny ; j++ ){
	x0 = w*(0.5+i)/nx;
	y0 = h*(0.5+j)/ny;
	
	selectImage(st_id);
	//Overlay.clear();
	//Overlay.drawString(""+i+" "+j,x0,y0);
	//Overlay.drawRect( x0-0.5*dx, y0-0.5*dy,dx,dy);
	//Overlay.show();	

	//Overlay.drawEllipse( x0-dx,y0-dy,2.*dx,2.*dy);
	//Overlay.show();

	// calculate the mode and width of the intensity distribution
	
	makeEllipse( x0-dx,y0-dy,x0+dx,y0+dy,1.);
	

	nmax=pow(2,bitDepth());
	print(nmax);
	if(nmax==256){
		getHistogram(vals,counts,256);
	}else{
		getHistogram(vals,counts,256,0,nmax);
	}
	tmp=get_hist_mode_lhm_hhm( counts );
	
	
	mode=vals[tmp[0]];
	lhm =mode-vals[tmp[1]];
	uhm =vals[tmp[2]]-mode;

	

	// copy the center tile 
	makeRectangle( x0-0.5*dx,y0-0.5*dy,dx,dy);

	run("Copy");
	run("Internal Clipboard");


	// Threshold
	run("16-bit");
	if( t_l < 0 ){
		setThreshold( 0,mode+t_l*lhm);
	}else{
		setThreshold( 0,mode+t_l*uhm);
	}
	//setThreshold( mode+t_l*uhm,160000);

	setOption("BlackBackground",false);
	run("Convert to Mask");
	run("Select All");
	run("Copy");
	close();


	selectImage(tt_id);
	makeRectangle(x0-0.5*dx, y0-0.5*dy,dx,dy);
	run("Paste");
	run("Invert");
	
}
}
selectImage(tt_id);
run("Set Scale...", "distance=1. known="+px2um+" unit=micron");
run("Make Binary");
//if( ! bmode_on){ setBatchMode(false);};
return(tt_id);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function measure_cells_in_image( id, ns ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
// Run the blob finder
selectImage(id);

Area_min = pi*size_min*size_min/4.;
Area_max = pi*size_max*size_max/4.;

run("Set Measurements...", "area centroid perimeter circularity fit bounding redirect=None decimal=3");

n0=nResults();
run("Analyze Particles...", "size="+Area_min+"-"+Area_max+"  circularity="+circ_min+"-"+circ_max+" record display add");
n1=nResults();

for(i=n0;i<n1;i++){
	setResult("Series",i,ns);
}

return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function display_counted_cells(n0,n1) {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
selectImage(w_id);
//setBatchMode(false);

paint_cells(n0,n1);
setLocation(Im_x0, Im_y0, Im_w0, Im_h0);
//setBatchMode(true);
return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function paint_cells(n0,n1){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Overlay.remove();
Overlay.clear();
setLineWidth(3);

setColor("yellow");
Overlay.drawRect(size_max/px2um,size_max/px2um,(width_um-2*size_max)/px2um,(height_um-2*size_max)/px2um);

if( n0 < n1 ){
for( j=n0; j<n1; j++){

	X=getResult("X",j)/px2um;
	Y=getResult("Y",j)/px2um;
	
	BX=getResult("BX",j)/px2um;
	BY=getResult("BY",j)/px2um;

	angle=(-getResult("Angle",j))*pi/180.;
	major=getResult("Major",j)/2./px2um;
	minor=getResult("Minor",j)/2./px2um;
	angle=getResult("Angle",j); 

	setColor(0, 0, 255);
	Overlay.drawString(""+(j+1),X,Y);
	exclude=getResult("Exclude",j);
	if(exclude==0){ // counted
		setColor("red");
		drawEllipse( X,Y,  major, minor, angle );
		Overlay.show();
	}
	
	if(exclude==1){	// excluded by size
	}	
	if(exclude==2){	// excluded by soft edge
		setColor("magenta");
		drawEllipse( X,Y,  major, minor, angle );
		Overlay.show();
	}		
	if(exclude==3){ // excluded by edge
		setColor("green");
		drawEllipse( X,Y,  major, minor, angle );
		Overlay.show();
	}			
}}

Overlay.show();
updateDisplay();

return;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function drawEllipse(x, y, a, b, angle) {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
  //    autoUpdate(false);
      setLineWidth(2);
      beta = -angle * (PI/180);
      for (i=0; i<=360; i+=12) {
          alpha = i*(PI/180) ;
          X = x + a*cos(alpha)*cos(beta) - b*sin(alpha)*sin(beta);
          Y = y + a*cos(alpha)*sin(beta) + b*sin(alpha)*cos(beta);
          if (i==0) Overlay.moveTo(X, Y); else Overlay.lineTo(X,Y);
          if (i==0)   {ax1=X; ay1=Y;}
          if (i==90)  {bx1=X; by1=Y;}
          if (i==180) {ax2=X; ay2=Y;}
          if (i==270) {bx2=X; by2=Y;}
      }
      //Overlay.drawLine(ax1, ay1, ax2, ay2);
      //Overlay.drawLine(bx1, by1, bx2, by2);
      updateDisplay;
}

// - - - - - - - - - - - - - - - - - - - - - - -
function get_meta_data( filepath ){
// - - - - - - - - - - - - - - - - - - - - - - -	
// Open the nd2 metadata and return the meta data string

path    =File.getParent(filepath);	     	// directory it is in
file    =substring(filepath,lengthOf(path)+1);

bio_options=" autoscale color_mode=Composite concatenate_series open_all_series view=Hyperstack stack_order=XYCZT";

// open the meta data
run("Bio-Formats Importer", "open=["+filepath+"] display_metadata view=[Metadata only]");
meta_win="Original Metadata - "+file;
selectWindow(meta_win);
meta_str=getInfo("window.contents");
run("Close");

return meta_str;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function extract_key_val( key, string ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	nl=lengthOf(string);
	nk=lengthOf(key);
	ns=indexOf( string, key);
	if( ns== -1 ){
		return "NA";
	}
	
	subs=substring(string,ns+nk);
	ne=indexOf(subs,"\n");
	val=substring(subs,0,ne);

	return val;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function test_and_close( name ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// checks is a window is open them closes it
if( isOpen(name) ){ 
	selectWindow(name); 
	run("Close"); 
}
return;
}


// - - - - - - - - - - - - - - - - - - - - - - - -
function list_files_in_dir( dir, ext){
// - - - - - - - - - - - - - - - - - - - - - - - -
// Returns a list of the files in a directory with the
// extension .ext
// if ext is "dir" it looks for directories
//
// Last updated: Wednesday, 20 August 2014 12:22:33

	// Check to see if we are looking for directories or files
	dir_mode=false;
	if( ext=="dir"){
		dir_mode=true;
	}

	// Get the files in the directoty
	files = getFileList(dir);

	// this is the internal extension
	if( lengthOf(ext)==0 || ext=="*" ){ // Wildcard matching
		ext_i="";
	}else{
		ext_i="."+ext;		
	}
	print("ext_i= "+ext_i);
	
	// Count how may of the files in appropriate extensions
	nfile=0;
	for (i=0; i<files.length; i++){
		file=files[i];

		if(dir_mode){
			if( File.isDirectory(dir+file)==1){
				nfile++;
			}
		}else{
			if( endsWith(file,ext_i) ){
				nfile++;
			}
		}
	}

	// create the holding array
	file_list=newArray(nfile);

	// go back through the list and copy the files into the array
	nf = 0;
	for (i=0; i<files.length; i++){
		file=files[i];
		if(dir_mode){
			if( File.isDirectory(dir+file)==1){
				file_list[nf]=file;
				nf++;
			}		
		}else{
			if( endsWith(file,ext_i) ){
				file_list[nf]=file;
				nf++;
			}
		}
	}
	return file_list;
}  


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function get_hist_mode_lhm_hhm( hist ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
// Input : 
// 	hist     histogram array
// Returns:
//	i : the bin containing the mode
nh=hist.length;

max_locs=Array.findMaxima(hist, 2);
mode = max_locs[0];	

hhalf=hist[mode]/2.;

chp=nh-1;
for(i=mode;i<nh;i++){
	if( hist[i]<hhalf ){ 
		chp=i;
		i=nh;
	}
}

chm=0;
for(i=mode;i>0;i--){
	if( hist[i]<hhalf ){ 
		chm=i;
		i=0;
	}
}

return newArray(mode,chm,chp);
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function get_dilution_factor_from_filename( string ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	df = 1.0;
	// remove the extension
	ftrim=trim_ext(file);
	
	// find instance of "d=", 
	ns=indexOf(ftrim,"d=");
	if( ns>1 ){
		// take stuff on right if it
		ftrim=substring(ftrim,ns+2);

		// trim the stuff to the right of it
		ne=indexOf(ftrim,"_");
		if(ne>0){
			ftrim=substring(ftrim,0,ne);
		}

		// parse the factor
		df=parseFloat(ftrim);
	}
	return df;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function trim_ext( string ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	nd=lastIndexOf(string,".");

	subs="";
	if(nd >=0 ){
		subs=substring(string,0,nd);
	}
	return subs; 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function get_ext( string ){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	nd=lastIndexOf(string,".");

	subs="";
	if(nd >=0 ){
		subs=substring(string,nd+1,lengthOf(string));
	}
	return subs; 
}



function erf(xi){
	x = xi;
    // save the sign of x
    sign=1;
    if( x < 1 ){ sign = -1;}
  
    x = abs(x);

    //constants
    a1 =  0.254829592;
    a2 = -0.284496736;
    a3 =  1.421413741;
    a4 = -1.453152027;
    a5 =  1.061405429;
    p  =  0.3275911;

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    return sign*y;
  } 
