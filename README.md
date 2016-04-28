# HySDeP_Hydrogen-Storage-Design-Platform

 Design and simulation platform for on-board hydrogen storage systems in Dymola, part of Andrea Mazzucco's PhD (Technical University of  Denmark, Oct. 2015).

 Installation: Make sure dymola is installed and working correctly (i.e., modelica and visual studio/express 2010 should be installed).

 Download all the files in "HySDeP_Hydrogen-Storage-Design-Platform" folder. Save the folder in Dymola folder in Documents. Extract the file contained in "ExternalMedia" into the "HySDeP_Hydrogen-Storage-Design-Platform" folder.

 Installation of CoolProp. Copy "CoolPropLib.lib" from CoolPropFiles folder into installation folder of dymola "Dymola/bin/lib". Copy  
 "CoolPropLib.h" from CoolPropFiles folder into installation folder of dymola "Dymola/source". Please note that the standard version of  CoolProp2Modelica does not work with this package, yet.

 Open Dymola, load "Hydrogen Storage Design Platform - HySDeP.mo" and then load "package.mo" (in: ExternalMedia\ExternalMedia 
 3.2.1) change directory to the base path of "HySDeP_Hydrogen-Storage-Design-Platform" and enjoy :) *


 *Note that in the component named "SAEJ6201" (in "CHG" component) you have to write the correct path for the lookupTables: APRR.mat, 
 FP.mat, SOC.mat and for HeatTransferProperties.txt
