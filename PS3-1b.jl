# ------------------------------------------------------------------- #
# Shu-Han Wang (sw2227)
# CHEME5440 Spring 2020
# Problem Set 3
# April 14, 2020
# ------------------------------------------------------------------- #

using LinearAlgebra
#Carbomoyl-p
stoichiometric_matrix[1,:] = [0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Citrulline
stoichiometric_matrix[2,:] = [-1 0 0 1 -2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Aspartate
stoichiometric_matrix[3,:] = [-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Argsucc.
stoichiometric_matrix[4,:] = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Fumarate
stoichiometric_matrix[5,:] = [0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0];
#Arginine
stoichiometric_matrix[6,:] = [0 1 -1 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Ornith.
stoichiometric_matrix[7,:] = [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
#Urea
stoichiometric_matrix[8,:] = [0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0];
#Ortho-phos.
stoichiometric_matrix[9,:] = [0 0 0 1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0];
#ATP
stoichiometric_matrix[10,:] = [-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
#AMP
stoichiometric_matrix[11,:] = [1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0];
#Diphos.
stoichiometric_matrix[12,:] = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0];
#H+
stoichiometric_matrix[13,:] = [0 0 0 0 3 -3 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
#O2
stoichiometric_matrix[14,:] = [0 0 0 0 4 -4 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0];
#NADPH
stoichiometric_matrix[15,:] = [0 0 0 0 3 -3 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
#NADP+
stoichiometric_matrix[16,:] = [0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0];
#NO
stoichiometric_matrix[17,:] = [0 0 0 0 -2 2 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0];
#H2O
stoichiometric_matrix[18,:] = [0 0 -1 0 -4 4 0 0 0 0 0 1 0 0 0 0 0 0 0 0 -1];


elemental_balance = Array{Float64,2}(undef,5,18)
elemental_balance[1,:]=[1 6 4 10 4 6 5 1 0 10 10 0 0 0 21 21 0 0]; #C
elemental_balance[2,:]=[4 13 7 18 4 14 12 4 3 16 14 4 1 0 30 29 0 2]; #H
elemental_balance[3,:]=[1 3 1 4 0 4 2 2 0 5 5 0 0 0 7 7 1 0]; #N
elemental_balance[4,:]=[5 3 4 6 4 2 2 1 4 13 7 7 0 2 17 17 1 1]; #O
elemental_balance[5,:]=[1 0 0 0 0 0 0 0 1 3 1 2 0 0 3 3 0 0]; #P


Balanced = elemental_balance*stoichiometric_matrix
# Results: First 6 columns = 0, thus stoichiometric_matrix is balanced!
