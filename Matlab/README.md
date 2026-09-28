# MATLAB Code

This is the code used to create and show the results in Royee Yosibash's
master thesis in electrical engineering: "IRREGULAR POLYNOMIAL CODES FOR
CODED COMPUTATION WITH NUMERICAL STABILITY AND GRACEFUL DEGRADATION".

All files tested and run in MATLAB version 2018a.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The project contains the following folders:

1. FramesTOOLBOX	- 	A custom "Toolbox" containing functions needed to run many of the scripts given in this project. It is ill-advised to change any of the
				functions in this Toolbox as it may lead to some scripts or functions in other folders crash or function unexpectedly. Most functions are
				either self explanitory or have documentation in the file, but some explanation are given to a few key functions:

				Shared code-construction, linear-algebra, eigenvalue, distribution, and
				statistics utilities belong in this folder. Workflow-specific classes
				and GUI-only helpers remain in their respective workflow folders.

				a. FrameParameters	-	This struct contains the parameters needed to define a frame type and dimensions. It currently supports a 
								multi-variable dimensions but only a single frame type in each instance construction. The struct also contains 
								multiple public functions that allow the user to create empirical statistics on the frame's numerical stability
								for multiple performance measures. While a constructor function exists for this struct, it is recomended using
								the GUI to create FameParameters objects.

				b. getCode 		-	This function creates a code generator matrix (if exists) that satisfies the parameters described by the inputs
								to the fucntion. Note that the codes implemented in this function were are only the codes that where necessity 
								during reasearch. Including new codes into this function is optional (and not difficult) but the reader should 
								by carefull, as changing the function carelessly might lead to unexpected results. 
								Note: The output of this function is the transpose ofthe generator matrix (this suits a frame-centric design
								conventions).

 
2. GUI 			-	This folder contains a GUI, that can be initiated using the "Run" command from the terminal or script. The app is designed in the MATLAB app 
				designer. The app is composed of three main areas:

				a. Frame creation area	-	The top of the app window allows the user to create frames of a desired code type. The user can choose to 
								create a frame with only one set of dimensions or hold on dimension constant and create multiple frames
								each with the other dimension taken from a grid (choosen by the user). The created frames apear in the table
								on the left hand side of the frame creation area. Note that "N" is the effectivly the number of rows of the 
								generator matrix and "M" the number of columns. 

				b. Functionality area	-	This area contains the many functions implemented in the app. The user can create plots using the selections in
								this area.

				c. Results plot area	-	The lower right hand side of the app is the graph in which the results are ploted into. The user can use the
								"copy axes" function in the Functionality area to copy the graph to an new external figure.

3. Coding Scheme	-	This folder contains the script and functions needed in order to compare coding schemes for the distributed computation of ther matrix-vector  
				and matrix-matrix cases. Follow the instructions on how to compare the codes under test and choose n,m,SNR and other parameters. The function
				also creates a TXT file that gives a table in LaTeX format.

## Setup

Before running an active workflow, initialize the project path by running
"setup" from the Matlab folder, or by calling "setup" after adding the Matlab
folder to the MATLAB path. This adds the active source folders without adding
the archived code. The active launch scripts call "setup" themselves and
resolve the project root from their own file locations.

## Active Entry Points

The current active entry points are:

1. GUI/Run.m			- Launches the Frame Analyzer GUI.
2. Coding Scheme/compareCodes.m	- Runs the coding-scheme comparison and writes result tables.
3. Results/CreateFiguresForArticle.m	- Creates the distribution comparison figures.

The default parameters for the coding-scheme comparison are defined in
Coding Scheme/compareCodesConfig.m. The function returns a configuration
struct and preserves the defaults used by compareCodes.m.





## TODO

1. Check that earlier versions of codes tested are all the transpose of the code generator matrix... 
2. Implement user friendly matrix-vector/matrix-matrix switch cases
