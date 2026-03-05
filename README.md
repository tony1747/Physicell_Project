#Add readme file to indicate the code changes in functions in cpp file.
Antigen threshold can be found in "threshold test.xlsx". Use it to update antigen threshold defined in PhysiCell_settings.xml based on corresponding samples.
Modify the introduce_il13_immune_cells() function in the cancer_immune_3D.cpp to update the injection loaction for IL13 CAR T-cells. Same applies to HER2 and EGFR CAR T-cells.
Modify the coloring function in the cancer_immune_3D.cpp to update the cell colorings for cancer cells, as they are different from samples to samples. "rgb color.txt" includes the RGB values for all three samples.
