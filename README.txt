Formation of synthetic RNA protein granules using engineered phage-coat-protein -RNA complexes

Supplamentray code prerequisits and instructions

The scripts and instructions provided were tested on Matlab 2021a 64-bit with the following additions:
- Curve Fitting Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

In addition, ImageJ FIJI (specifically, the Mosaic plugin) was used to read microscopy image files, detect bright spots, and convert to the formats required for the scripts.
All scripts were built and tested in a Microsoft Windows 10 version 21H1 environment.

Installation instructions: Unzip the file and open the scripts in Matlab
Expected running time of all the scripts: minutes on an average desktop computer

Microscopy image files:
All microscopy files analyzed using the provided scripts were captured on a Nikon Eclipse Ti-E epifluorescent microscope with a 100x 1.45 NA oil immersion objective.
Images were captured every 10 seconds for a duration of 1 hour. 
Sample files provided as a demo in the folder '8x' depicting PCP-8x granules tracked in two channles (files were converted from the Nikon proprietary ND2 to tiff using imageJ)

Provided Matlab scripts (indentation means sub-function used by a script)
- read_tables.m
- collect_intensity_data.m
	- generate_movArray.m
- collect_statistical_data.m
	- get_position.m
		- runLengthEncode.m
	- plot_tracks_custom.m
- plot_2_channel_measurements.m
	- plot_tracks_custom_subplot_version.m
- fit_data_to_poisson.m
	- poisson_sub_func.m
	- violin.m
- Generate_plots_for_paper.m


The code is divided into 4 standalone parts:

- Full analysis pipeline (scripts read_tables, collect_intensity_data, collect_statistical_data)
	Intensity signal gathering and identification of molecule entry and exit events.
	Detailed instructions provided in 'Full pipeline walkthrough.pdf'
- Fitting to Poisson (script fit_data_to_poisson)
	Fitting of intensity data (measurements of median intensity, or burst amplitude measurements) to a modified Poisson distribution
	Detailed isntructions provided in 'Plot fits to Poisson.pdf'
- Plot measurements from 2 channels
	Ploting of intensity signals from granules measured in two fluorescence channles (e.g., 488 [nm] and 585 [nm]). 
	Detailed instructions provided in 'Plot 2 channel measurements.pdf'
- Recreate plots from the paper
	Plotting of the figures appearing in the paper
	Script can be run as provided 

