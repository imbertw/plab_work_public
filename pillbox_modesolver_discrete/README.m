%README: This doc explains what each of the codes are doing in this
%directory
%% Analytical Pillbox Modesolver
%1. CylindricalParameter: Calculates the propagation constant for a
%cylindrical waveguide
%2. PillboxFields: Given the values from Cylindrical Parameter, computes
%the fields of the pillbox
%3. PillboxIndex: Calculates what the index of the imaginary index material
%in the corners should be
%4. PillboxLT: Calculates the length of the pillbox mode, and computes the
%longitudinal field E(z).
%5. PillboxTV: Calculates the propagation constants and computes the
%transverse field E(rho) of the pillbox mode. Uses CylindricalParameter and
%PillboxFields.
%6. PillboxMS: The pillbox modesolver. Uses PillboxIndex, PillboxLT,
%PillboxTV, PillboxMV, ThreeDField, ThreeDFieldPlotIsosurface to generate
%the plots of the pillbox mode.
%7. PillboxMV: Calculates the mode volume of the Pillbox mode
%8. ThreeDField: Computes the 3D Field from inputs from PillboxTV and PillboxLT 
%9. ThreeDFieldPlotIsosurface: Generates an Isosurface plot of the 3D
%field. Might be useful for figures in paper.

%% Discretized Pillbox Modesolver
%1. s_m2cylphiR1: Calculates propagation constant and computes E(rho) for a
%cylindrical waveguide, and compares it to the pillbox mode from the
%analytical modesolver.
%2. cylindrical_modesolver.m: Imbert's first attempt at writing a
%cylindrical modesolver, from 2017, during Ramachandran's class on fiber
%waveguides.