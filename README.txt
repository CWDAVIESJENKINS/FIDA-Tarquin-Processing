CWDaviesJenkins
Set of functions for processing of Tarquin (Wilson étal 2011) data via Matlab. Written at CUBRIC, Cardiff University under WSA grant. Uses FID-A  (https://github.com/CIC-methods/FID-A), and structural functions were adapted from co-registration functions in Gannet (https://github.com/richardedden/Gannet3.1). A modified version of FID-A may be found in lib/

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

