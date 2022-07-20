# define gravity and pi 
set g 9.80665 
set pi [expr atan(1)*4] 

#---SOIL GEOMETRY 

# thicknesses of soil profile (m) 
set soilThick 20.21 

# number of soil layers (not counting layer 0 which is bedrock) 
set numLayers 9 

# layer thicknesses 
set layerThick(9) 1.4 
set layerThick(8) 2.1 
set layerThick(7) 1.5 
set layerThick(6) 1.3 
set layerThick(5) 1.9 
set layerThick(4) 1.275 
set layerThick(3) 1.275 
set layerThick(2) 5.02 
set layerThick(1) 3.44 
set layerThick(0) 1.0 

# depth of water table (create a new layer at WT) 
# if water not present set waterTable anywhere below depth of model 
set waterTable 3.5 

#---BASIC MATERIAL PROPERTIES 

# soil mass density (Mg/m^3) 
set rho(9) 1.9235474006116209 
set rho(8) 2.0305810397553516 
set rho(7) 1.9347604485219163 
set rho(6) 2.058103975535168 
set rho(5) 1.9490316004077473 
set rho(4) 1.8807339449541283 
set rho(3) 1.8430173292558611 
set rho(2) 1.8348623853211008 
set rho(1) 1.9164118246687054 
set rho(0) 2.038735983690112 
set rhoWater 1.0 

# soil shear wave velocity for each layer (m/s) 
set Vs(9) 110.0 
set Vs(8) 195.0 
set Vs(7) 195.0 
set Vs(6) 195.0 
set Vs(5) 175.0 
set Vs(4) 165.0 
set Vs(3) 165.0 
set Vs(2) 197.0 
set Vs(1) 273.0 
set Vs(0) 500.0 


# Poisson ratio of soil 
set nu 0.25 

# rock elastic properties) 
# bedrock shear wave velocity (m/s) 
set rockVS $Vs(0) 
# bedrock mass density (Mg/m^3) 
set rockDen $rho(0) 


#---OTHER PDMY02 MODEL PARAMETERS 

# soil friction angle (°) 
set phi(9) 38.0 
set phi(8) 40.0 
set phi(7) 37.0 
set phi(6) 40.0 
set phi(5) 40.0 
set phi(4) 35.0 
set phi(3) 35.0 
set phi(2) 35.0 
set phi(1) 37.0 
set phi(0) 40.0 

# peak shear strain 
set gammaPeak 0.1 

# relative density (%) 
# Layer 9: Dr = 71.0 
# Layer 8: Dr = 99.5 
# Layer 7: Dr = 89.6 
# Layer 6: Dr = 99.4 
# Layer 5: Dr = 59.2 
# Layer 4: Dr = 29.7 
# Layer 3: Dr = 19.85 
# Layer 2: Dr = 40.0 
# Layer 1: Dr = 65.0 
# Layer 0: Dr = 100.0 

# reference pressure is computed as mean confining pressure 
# at refDepth for each layer (0 is ToL, 1 is BoL) 
set refDepth(9) 0.5 
set refDepth(8) 1.0 
set refDepth(7) 0.5 
set refDepth(6) 0.5 
set refDepth(5) 0.5 
set refDepth(4) 0.5 
set refDepth(3) 0.0 
set refDepth(2) 0.5 
set refDepth(1) 0.5 
set refDepth(0) 0.5 

# pressure dependency coefficient 
set pressCoeff(9) 0.25 
set pressCoeff(8) 0.25 
set pressCoeff(7) 0.0 
set pressCoeff(6) 0.0 
set pressCoeff(5) 0.0 
set pressCoeff(4) 0.0 
set pressCoeff(3) 0.5 
set pressCoeff(2) 0.5 
set pressCoeff(1) 0.5 
set pressCoeff(0) 0.0 

# phase transformation angle (not for layer 0) (°) 
set phaseAng(9) 26.0681 
set phaseAng(8) 26.5 
set phaseAng(7) 26.49032 
set phaseAng(6) 26.5 
set phaseAng(5) 26.20748 
set phaseAng(4) 31.0 
set phaseAng(3) 31.0 
set phaseAng(2) 30.0 
set phaseAng(1) 26.3 

# contraction (not for layer 0) 
set contract1(9) 0.0194 
set contract1(8) 0.016 
set contract1(7) 0.01568 
set contract1(6) 0.016 
set contract1(5) 0.04264 
set contract1(4) 0.087 
set contract1(3) 0.087 
set contract1(2) 0.067 
set contract1(1) 0.032 
set contract2(9) 1.4931 
set contract2(8) 1.4 
set contract2(7) 1.45032 
set contract2(6) 1.4 
set contract2(5) 3.4383199999999996 
set contract2(4) 5.0 
set contract2(3) 5.0 
set contract2(2) 4.5 
set contract2(1) 2.1 
set contract3(9) 0.1485 
set contract3(8) 0.14 
set contract3(7) 0.1392 
set contract3(6) 0.14 
set contract3(5) 0.21503999999999998 
set contract3(4) 0.3 
set contract3(3) 0.3 
set contract3(2) 0.27 
set contract3(1) 0.18000000000000002 

# dilation coefficients (not for layer 0) 
set dilate1(9) 0.16349999999999998 
set dilate1(8) 0.25 
set dilate1(7) 0.24719999999999998 
set dilate1(6) 0.25 
set dilate1(5) 0.06736 
set dilate1(4) 0.01 
set dilate1(3) 0.01 
set dilate1(2) 0.02 
set dilate1(1) 0.10200000000000001 
set dilate2(9) 3.0 
set dilate2(8) 3.0 
set dilate2(7) 3.0 
set dilate2(6) 3.0 
set dilate2(5) 3.0 
set dilate2(4) 3.0 
set dilate2(3) 3.0 
set dilate2(2) 3.0 
set dilate2(1) 3.0 
set dilate3(9) 0.0 
set dilate3(8) 0.0 
set dilate3(7) 0.0 
set dilate3(6) 0.0 
set dilate3(5) 0.0 
set dilate3(4) 0.0 
set dilate3(3) 0.0 
set dilate3(2) 0.0 
set dilate3(1) 0.0 

# void ratio (need it for layer 0 for element definition) 
set voidR(9) 0.6404 
set voidR(8) 0.58 
set voidR(7) 0.5808800000000001 
set voidR(6) 0.58 
set voidR(5) 0.67148 
set voidR(4) 0.76 
set voidR(3) 0.76 
set voidR(2) 0.73 
set voidR(1) 0.656 
set voidR(0) 0.58 

