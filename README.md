# HySDeP_Hydrogen-Storage-Design-Platform

# Design and simulation platform for on-board hydrogen storage systems in Dymola, part of Andrea Mazzucco's PhD (Technical University of # Denmark, Oct. 2015).

# Installation: Make sure dymola is installed and working correctly (i.e., modelica and visual studio/express 2010 should be installed).

# Download HySDeP package. Extract files in Dymola folder in Documents. Keep library structure within the "HySDeP" folder.

# Installation of CoolProp. Copy "CoolPropLib.lib" from CoolPropFiles folder into installation folder of dymola "Dymola/bin/lib". Copy  
# "CoolPropLib.h" from CoolPropFiles folder into installation folder of dymola "Dymola/source". Please note that the standard version of # CoolProp2Modelica does not work with this package, yet.

# Open Dymola, load "Hydrogen Storage Design Platform - HySDeP.mo" and then load "package.mo" (in: ExternalMedia\Modelica\ExternalMedia 
# 3.2.1) change directory to the base path of "HySDeP" and enjoy :) *


# *Note that in the component named "SAEJ6201" (in "CHG" component) you have to write the correct path for the lookupTables: APRR.mat, 
# FP.mat, SOC.mat and HeatTransferProperties.txt
