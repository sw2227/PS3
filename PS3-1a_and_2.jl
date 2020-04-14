# ------------------------------------------------------------------- #
# Shu-Han Wang (sw2227)
# CHEME5440 Spring 2020
# Problem Set 3
# April 14, 2020
# ------------------------------------------------------------------- #

include("Flux.jl")

using LinearAlgebra
metabolite_concentration = Array{Float64,2}(undef,18,1)
metabolite_concentration[1] = 0.0; #Carbomoyl-p
metabolite_concentration[2] = 0.0; #Citrulline
metabolite_concentration[3] = 0.0149; #Aspartate
metabolite_concentration[4] = 0.0; #Argsucc.
metabolite_concentration[5] = 0.0004850; #Fumarate
metabolite_concentration[6] = 0.000255; #Arginine
metabolite_concentration[7] = 0.004489; #Ornith.
metabolite_concentration[8] = 0.0; #Urea
metabolite_concentration[9] = 0.0; #Ortho-phos.
metabolite_concentration[10] = 0.00467; #ATP
metabolite_concentration[11] = 0.0000423; #AMP
metabolite_concentration[12] = 0.0; #Diphos.
metabolite_concentration[13] = 0.0; #H+
metabolite_concentration[14] = 0.0; #O2
metabolite_concentration[15] = 0.0000654; #NADPH
metabolite_concentration[16] = 0.0000284; #NADP+
metabolite_concentration[17] = 0.0; #NO
metabolite_concentration[18] = 0.0; #H2O
metabolite_uconc = metabolite_concentration.*1000000
cell_mass = 2.3e-9;
cell_volume = 1.0e-12;
water_X = 0.798;
celldrymass = (1 - water_X)*cell_mass;
convert(uM) = uM*(cell_volume)/celldrymass;
metabolite = convert(metabolite_uconc)


stoichiometric_matrix = Array{Float64,2}(undef,18,21)
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


kcat = Array{Float64,2}(undef,6,1)
kcat[1] = 203.0;
kcat[2] = 34.5;
kcat[3] = 249.0;
kcat[4] = 88.1;
kcat[5] = 13.7;
kcat[6] = 13.7;


kM = Array{Float64,2}(undef,11,1)
kM[1] = 0.00013; #Carbomoyl-p
kM[2] = 0.000056; #Citrulline1
kM[3] = 0.0; #Citrulline2
kM[4] = 0.000154; #Aspartate
kM[5] = 0.0001; #Argsucc.
kM[6] = 0.00155; #Arginine1
kM[7] = 3.497; #Arginine2
kM[8] = 0.0016; #Ornith.
kM[9] = 0.000392; #ATP
kM[10] = 0.0000003; #NADPH
kM[11] = 0.0; #NADP+
ukM = kM.*1000000;
kMa = convert(ukM)


E = 0.01;
seta = 1;
#Lower Bounds
default_bounds_array = zeros(21,2);
UNIT = 10000/3600
default_bounds_array[11,1] = -(UNIT);
default_bounds_array[12,1] = -(UNIT);
default_bounds_array[13,1] = -(UNIT);
default_bounds_array[14,1] = -(UNIT);
default_bounds_array[15,1] = -(UNIT);
default_bounds_array[16,1] = -(UNIT);
default_bounds_array[17,1] = -(UNIT);
default_bounds_array[18,1] = -(UNIT);
default_bounds_array[19,1] = -(UNIT);
default_bounds_array[20,1] = -(UNIT);
default_bounds_array[21,1] = -(UNIT);
#Upper Bounds
default_bounds_array[1,2] = E*seta*kcat[1]*
(metabolite[3]/(metabolite[3]+kMa[4]))*(metabolite[10]/(metabolite[10]+kMa[9]));
default_bounds_array[2,2] = E*seta*kcat[2];
default_bounds_array[3,2] = E*seta*kcat[3]*
(metabolite[6]/(metabolite[6]+kMa[6]));
default_bounds_array[4,2] = E*seta*kcat[4]*
(metabolite[7]/(metabolite[7]+kMa[8]));
default_bounds_array[5,2] = E*seta*kcat[5];
default_bounds_array[6,2] = E*seta*kcat[6]*
(metabolite[6]/(metabolite[6]+kMa[7]))*
(metabolite[15]/(metabolite[15]+kMa[10]));
default_bounds_array[7,2] = (UNIT);
default_bounds_array[8,2] = (UNIT);
default_bounds_array[9,2] = (UNIT);
default_bounds_array[10,2] = (UNIT);
default_bounds_array[11,2] = (UNIT);
default_bounds_array[12,2] = (UNIT);
default_bounds_array[13,2] = (UNIT);
default_bounds_array[14,2] = (UNIT);
default_bounds_array[15,2] = (UNIT);
default_bounds_array[16,2] = (UNIT);
default_bounds_array[17,2] = (UNIT);
default_bounds_array[18,2] = (UNIT);
default_bounds_array[19,2] = (UNIT);
default_bounds_array[20,2] = (UNIT);
default_bounds_array[21,2] = (UNIT);


species_bounds_array = zeros(18,2);
objective_coefficient_array = zeros(21);
objective_coefficient_array[10] = -1;
OPTFLUX = calculate_optimal_flux_distribution(stoichiometric_matrix,
default_bounds_array,species_bounds_array,objective_coefficient_array);
Z = 3600*OPTFLUX[2]/1000;
urea_only = Z[10];

print("-------------------------------------------------------------------\n")
print("The maximum rate of urea production is $(urea_only) mmol/gDW-hr. \n")
print("------------------------------------------------------------------- \n")
print("Stoichiometric matrix: \n")
println(stoichiometric_matrix)
