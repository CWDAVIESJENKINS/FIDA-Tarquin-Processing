CWDaviesJenkins
Set of functions for processing of Tarquin (Wilson étal 2011) data via Matlab. Written at CUBRIC, Cardiff University under WSA grant. Uses a modified FID-A (https://github.com/CIC-methods/FID-A), and structural functions were adapted from co-registration functions in Gannet along with some other functions (https://github.com/richardedden/Gannet3.1).

Before running, please update Tarquin_Config.m to reflect your Tarquin installation path.


Quick summary of available functions:

Configuration of options:
	Tarquin_Config.m

Main run function:
	Tarquin_Run

Read Tarquin csv data into Matlab
	Tarquin_Read.m

Process stats of cell array of Tarquin fit results:
	Tarquin_Process.m	

Plotting tools:
	Tarquin_Plot.m
	Tarquin_Plot2.m
	Tarquin_Checker.m

Functions based on Gannet tissue correction methods:
	Struct_CoRegister.m
	Struct_Segment.m
	Struct_Quantify.m	
	Struct_AddSegment.m	
	Struct_VoxOverlap.m

Attempt to generate parameter file from FID-A format
	Tarquin_GenPara.m	

Read TWIX for Tarquin
	Tarquin_ReadDat.m
	Tarquin_ReadDat_RSR.m

Please cite the following:

Wilson, Martin, et al. "A constrained least‐squares approach to the automated quantitation of in vivo 1H magnetic resonance spectroscopy data." Magnetic resonance in medicine 65.1 (2011): 1-12.

Simpson, Robin, et al. "Advanced processing and simulation of MRS data using the FID appliance (FID‐A)—an open source, MATLAB‐based toolkit." Magnetic resonance in medicine 77.1 (2017): 23-33.

Edden, Richard AE, et al. "Gannet: A batch‐processing tool for the quantitative analysis of gamma‐aminobutyric acid–edited MR spectroscopy spectra." Journal of Magnetic Resonance Imaging 40.6 (2014): 1445-1452.

Also, if using Violin plots:
Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany. hhoffmann@uni-bonn.de