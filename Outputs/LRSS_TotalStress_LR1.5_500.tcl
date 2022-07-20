########################################################### 
#                                                         # 
# Total Stress Site Response Analysis                     # 
# Constitutive Models: PDMY02 (sands) & PIMY (clays)      # 
# Pore Pressure Generation Allowed: No                    # 
#                                                         # 
########################################################### 

# Extract exterior inputs and give them a variable name 
set site [lindex $argv 0] 
set modelID [lindex $argv 1] 
set gMotionPath [lindex $argv 2] 
set gMotionName [lindex $argv 3] 
set npts [lindex $argv 4] 
set dt [lindex $argv 5] 
set saveDir [lindex $argv 6] 

set dash "_" 
set analID $site$dash$modelID$dash$gMotionName
set siteModelID $site$dash$modelID

# Define gravity and pi 
set g 9.80665 
set pi [expr atan(1)*4] 

# ----------------------------------------------------------------------------------------- 
#  1. DEFINE SOIL GEOMETERY AND MATERIAL PARAMETERS 
# ----------------------------------------------------------------------------------------- 

#---SOIL GEOMETRY 

# Thicknesses of soil profile (m) 
set soilThick 20.21 

# Number of soil layers (not counting layer 0 which is bedrock) 
set numLayers 9 

# Layer thicknesses 
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

# Depth of water table (create a new layer at WT) 
# If water not present set waterTable anywhere below depth of model 
set waterTable 3.5 

#---BASIC MATERIAL PROPERTIES 

# Soil mass density (Mg/m^3) 
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

# Soil shear wave velocity for each layer (m/s) 
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

# Rock elastic properties) 
# Bedrock shear wave velocity (m/s) 
set rockVS $Vs(0) 
# Bedrock mass density (Mg/m^3) 
set rockDen $rho(0) 

#--- MODEL PARAMETERS 

# Consitutive model 
# constModelFlag: 1 for PDMY02 (sands), 0 for PIDMY (clays) 
set constModelFlag(9) 1 
set constModelFlag(8) 1 
set constModelFlag(7) 1 
set constModelFlag(6) 1 
set constModelFlag(5) 1 
set constModelFlag(4) 1 
set constModelFlag(3) 1 
set constModelFlag(2) 1 
set constModelFlag(1) 1 
set constModelFlag(0) 1 

# Soil friction angle (°) 
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

# Peak shear strain 
set gammaPeak 0.1 

# Relative density (%) 
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

# Reference pressure is computed as mean confining pressure 
# at refDepth for each layer (0 is ToL, 1 is BoL, 0.5 is in the middle) 
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

# Pressure dependency coefficient 
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

# Phase transformation angle (not for layer 0) (°) 
set phaseAng(9) 26.0681 
set phaseAng(8) 26.5 
set phaseAng(7) 26.49032 
set phaseAng(6) 26.5 
set phaseAng(5) 26.20748 
set phaseAng(4) 31.0 
set phaseAng(3) 31.0 
set phaseAng(2) 30.0 
set phaseAng(1) 26.3 

# Contraction coefficients (not for layer 0) 
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

# Dilation coefficients (not for layer 0) 
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

# Void ratio (need it for layer 0 for element definition) 
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

# Constitutive model flag. 1 for PDMY02 (sands), 0 for PIMY (clays)
array set constModelFlag [list 4 1 3 1 2 1 1 1] 

# Allow excess pore pressure generation? Yes or No
# If No, permeability is automatically set very high for dynamic analysis
set allowPWP       No

# Define layer boundaries
set layerBound(0) $layerThick(0)
puts "layer boundary 0 = $layerBound(0)"
for {set i 1} {$i <= $numLayers} {incr i 1} {
    set layerBound($i) [expr $layerBound([expr $i-1])+$layerThick($i)]
    puts "layer boundary $i = $layerBound($i)"
}

# Soil cohesion for clay layers (i.e., PIMY layers)        
# array set cohesion [list 4 40.0 2 65.0]

# Flags for water table and Vs inversion
# set VsInvTopLayer to "Yes" if there is a velocity inversion immediately below upper layer (else "No")
# set waterTopLayer to "Yes" if water table within upper most layer and layer was split in two (else "No")
# if waterTopLayer == "Yes", should set refDepth(numLayers) = 1.0 and refDepth(numLayers-1) = 0.0
set VsInvTopLayer "No"
set waterTopLayer "No"

# Computed as mean confining pressure at refDepth for each layer (0 is ToL, 1 is BoL)

## for top layer
if {$layerThick($numLayers) > $waterTable} {
    set vertStress($numLayers) [expr ($rho($numLayers) - $rhoWater) * $g * $layerThick($numLayers) * $refDepth($numLayers)]
}   else {
    set vertStress($numLayers) [expr $rho($numLayers) * $g * $layerThick($numLayers) * $refDepth($numLayers)]
}
set kNot($numLayers) [expr (1.0 - sin($phi($numLayers) * 2.0 * $pi / 360.0))]
set meanStress [expr $vertStress($numLayers) * (1.0 + 2.0 * $kNot($numLayers)) / 3.0]
set refPress($numLayers) $meanStress

## for other layers (not necessary for layer 0 bc it's lin elastic bedrock)
set bottomDepth($numLayers) $layerThick($numLayers)

for {set k [expr $numLayers - 1]} {$k > 0 && $k <= [expr $numLayers - 1]} {incr k -1} {
    set bottomDepth($k) [expr $bottomDepth([expr $k + 1]) + $layerThick($k)]
    if {$bottomDepth($k) > $waterTable && $bottomDepth([expr $k + 1]) > $waterTable} {
        set vertStress($k) [expr $vertStress([expr $k + 1]) + ($rho([expr $k + 1]) - $rhoWater) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + ($rho($k) - $rhoWater) * $g * $layerThick($k) * $refDepth($k)]
} elseif {$bottomDepth($k) > $waterTable} {
    set vertStress($k) [expr $vertStress([expr $k + 1]) + $rho([expr $k + 1]) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + ($rho($k) - $rhoWater) * $g * $layerThick($k) * $refDepth($k)]
} else {
    set vertStress($k) [expr $vertStress([expr $k + 1]) + $rho([expr $k + 1]) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + $rho($k) * $g * $layerThick($k) * $refDepth($k)]
}
set kNot($k) [expr 1.0 - sin($phi($k) * 2.0 * $pi / 360.0)]
set meanStress [expr $vertStress($k) * (1.0 + 2.0 * $kNot($k)) / 3.0]
set refPress($k) $meanStress
}

# Compute Vs_not, the constant required for an exponential function of Vs to have
# equal travel time to a layer of constant Vs
set Vs_not [expr $Vs($numLayers)*pow($layerThick($numLayers), -$pressCoeff($numLayers) / 2.0) / (1.0 - $pressCoeff($numLayers) / 2.0)]
puts "Vs_not = $Vs_not"

# Soil shear modulus for each layer (kPa)

# for top layer
if {$VsInvTopLayer == "Yes" || $waterTopLayer == "Yes"} {
    set G($numLayers) [expr $rho($numLayers)*$Vs($numLayers)*$Vs($numLayers)]
} else {
    set G($numLayers) [expr $rho($numLayers) * pow($Vs_not, 2.0 ) * pow($rho($numLayers) * $g  * (1.0 + 2.0 * $kNot($numLayers)) / (3.0 * $refPress($numLayers)), -$pressCoeff($numLayers)) ]
}
# for all other layers    
for {set k 0} {$k < $numLayers} {incr k 1} {
    set G($k) [expr $rho($k)*$Vs($k)*$Vs($k)]    
}

# Soil elastic modulus for each layer (kPa)
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set E($k)       [expr 2*$G($k)*(1+$nu)]
}
# Soil bulk modulus for each layer (kPa)
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set bulk($k)    [expr $E($k)/(3*(1-2*$nu))]
}

# Print the value of Shear Modulus to output
if {[file exists "$saveDir/$siteModelID.nodeInfo.txt"] == "0"} {
    parray G
} 

#-----------------------------------------------------------------------------------------
#  2. MESH GEOMETRY
#-----------------------------------------------------------------------------------------

#---MESH GEOMETRY
# highest frequency desired to be well resolved (Hz)
set fMax    25.0
# Wavelength of highest resolved frequency for each layer
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set wave($k) [expr $Vs($k)/$fMax]
}
# Number of elements per one wavelength
set nEle     8

# Determine number of elements in column
set nElemT 0
for {set k 0} {$k <= $numLayers} {incr k 1} {

    # maximum element height
    set hEleMax($k) [expr $wave($k)/$nEle]

    # interger number of elements
    set nElemY($k) [expr int(floor($layerThick($k)/$hEleMax($k))+1)]
    puts "number of vertical elements in layer $k: $nElemY($k)"

    set nElemT [expr $nElemT + $nElemY($k)]

    # actual element height
    set sElemY($k) [expr {$layerThick($k)/$nElemY($k)}] 

    puts "vertical size of elements in layer $k: $sElemY($k)"
}
puts "total number of vertical elements: $nElemT"

# Number of nodes in vertical direction
set nNodeY  [expr 2*$nElemT+1]

# Number of elements in horizontal direction
set nElemX  1
# Number of nodes in horizontal direction
set nNodeX  [expr 2*$nElemX+1]
# Horizontal element size (m)
set sElemX  2.0

# Total number of nodes
set nNodeT  [expr $nNodeX*$nNodeY]

#-----------------------------------------------------------------------------------------
#  3. CREATE PORE PRESSURE NODES AND FIXITIES
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 2 -ndf 3

set ppNodesInfo [open $saveDir/$siteModelID.ppNodesInfo.dat w]
set LeftPPNodesInfo [open $saveDir/$siteModelID.LeftPPNodesInfo.dat w]
set dryNodeCount 1
set PPNodeCount 1
set layerNodeCount 0
set yCoordCount 0
# loop over soil layers
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in horizontal direction
    for {set i 1} {$i <= $nNodeX} {incr i 2} {
        # loop in vertical direction
        if {$k == 0} {
            set bump 1
    } else {
        set bump 0
    }
    for {set j 1} {$j <= [expr 2*$nElemY($k)+$bump]} {incr j 2} {

        set xCoord       [expr ($i-1)*$sElemX/2.0]
        set yctr    [expr $j + $layerNodeCount]
        if {$k == 0} {
            set yCoord       [expr ($j-1)*$sElemY($k)/2.0]
      } else {
          set yCoord       [expr $layerBound([expr $k - 1]) + $sElemY($k) + ($j - 1)*$sElemY($k)/2.0]
      }
      set nodeNum      [expr $i + ($yctr-1)*$nNodeX]

      node $nodeNum  $xCoord  $yCoord

      # puts "yctr = $yctr"
      # puts "xCoord = $xCoord"
      # puts "yCoord = $yCoord"
      # puts "nodeNum = $nodeNum"

      set PPNode($PPNodeCount) $nodeNum
      set PPNodeCount [expr $PPNodeCount + 1]

      # output nodal information to data file
      puts $ppNodesInfo "$nodeNum  $xCoord  $yCoord"

      if {$xCoord == 0} {
          puts $LeftPPNodesInfo "$nodeNum  $xCoord  $yCoord"
          #puts "this is a PP Node"
      }

      # designate nodes above water table
      set waterHeight [expr $soilThick-$waterTable]
      if {$yCoord>=$waterHeight} {
          set dryNode($dryNodeCount) $nodeNum
          set dryNodeCount [expr $dryNodeCount+1]
        }
    }
}
set layerNodeCount [expr $yctr + 1]
}
close $ppNodesInfo
close $LeftPPNodesInfo
puts "Finished creating all -ndf 3 nodes..."

# define fixities for pore pressure nodes above water table
for {set i 1} {$i < $dryNodeCount} {incr i 1} {
    fix $dryNode($i)  0 0 1
}

# define fixities for pore pressure nodes at base of soil column
fix 1  0 1 0
fix 3  0 1 0
puts "Finished creating all -ndf 3 boundary conditions..."


# define equal degrees of freedom for pore pressure nodes
for {set i 1} {$i <= [expr 3*$nNodeY-2]} {incr i 6} {
    equalDOF $i [expr $i+2]  1 2
}
puts "Finished creating equalDOF for pore pressure nodes..."

if {[file exists "$saveDir/$siteModelID.PPnodeInfo.txt"] == "0"} {
    print $saveDir/$siteModelID.PPnodeInfo.txt -node
} 
#-----------------------------------------------------------------------------------------
#  4. CREATE INTERIOR NODES AND FIXITIES
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 2 -ndf 2

# central column of nodes
set xCoord  [expr $sElemX/2]
# loop over soil layers
set layerNodeCount 0
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in vertical direction
    if {$k == 0} {
        set bump 1
} else {
    set bump 0
}
for {set j 1} {$j <= [expr 2 * $nElemY($k) + $bump]} {incr j 1} {

    set yctr    [expr $j + $layerNodeCount]
    if {$k == 0} {
        set yCoord  [expr ($j - 1) * $sElemY($k) / 2.0]
  } else {
      set yCoord  [expr $layerBound([expr $k - 1]) + $sElemY($k) / 2.0 +  ($j - 1) * $sElemY($k) / 2.0]
  }      
  set nodeNum [expr 3 * $yctr - 1] 

  node  $nodeNum  $xCoord  $yCoord 
}
set layerNodeCount $yctr
}

# Interior nodes on the element edges
# loop over layers
set layerNodeCount 0
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in vertical direction
    for {set j 1} {$j <= $nElemY($k)} {incr j 1} {

        set yctr [expr $j + $layerNodeCount]
        if {$k == 0} {
            set yCoord   [expr $sElemY($k) * ($j - 0.5)]
  } else {
      set yCoord   [expr $layerBound([expr $k - 1]) + ($j - 0.5) * $sElemY($k)]
  }
  set nodeNumL [expr 6*$yctr - 2]
  set nodeNumR [expr $nodeNumL + 2]

  node  $nodeNumL  0.0  $yCoord
  node  $nodeNumR  $sElemX  $yCoord
}
set layerNodeCount $yctr
}
puts "Finished creating all -ndf 2 nodes..."

# Define fixities for interior node at base of soil column
fix 2  0 1
puts "Finished creating all -ndf 2 boundary conditions..."

# Define equal degrees of freedom which have not yet been defined
for {set i 1} {$i <= [expr 3*$nNodeY-6]} {incr i 6} {
    equalDOF $i          [expr $i+1]  1 2
    equalDOF [expr $i+3] [expr $i+4]  1 2
    equalDOF [expr $i+3] [expr $i+5]  1 2
}
equalDOF [expr $nNodeT-2] [expr $nNodeT-1]  1 2
puts "Finished creating equalDOF constraints..."

# Print all node info to text file (if statement so it only does it once when batching)
puts "ground motion is $gMotionName"
if {[file exists "$saveDir/$siteModelID.nodeInfo.txt"] == "0"} {
    print $saveDir/$siteModelID.nodeInfo.txt -node
} 

#-----------------------------------------------------------------------------------------
#  5. CREATE SOIL MATERIALS
#-----------------------------------------------------------------------------------------

# Define grade of slope (%)
set grade 0.0
set slope [expr atan($grade/100.0)]
set g -9.81
set bulkWater 2.2e6
set materialProperties [open $saveDir/$siteModelID.materialProperties.dat w]

# Define nonlinear material for soil
for {set i 1} {$i <= $numLayers} {incr i 1} {

    if {$constModelFlag($i) == 1} {
	
		# PDMY02 Model
        nDMaterial PressureDependMultiYield02 $i 2 $rho($i) $G($i) $bulk($i) $phi($i) $gammaPeak \
        $refPress($i) $pressCoeff($i) $phaseAng($i) \
        $contract1($i) $contract3($i) $dilate1($i) $dilate3($i) \
        20 $contract2($i) $dilate2($i) 1 0 $voidR($i) 0.9 0.02 0.7 101.0  

        set thick($i) 1.0
        set xWgt($i)  [expr $g*sin($slope)]
        set yWgt($i)  [expr $g*cos($slope)]
        set porosity($i) [expr $voidR($i) / (1 + $voidR($i))]
        set uBulk($i) [expr $bulkWater/$porosity($i)]
        set hPerm($i) 1.0e-4
        set vPerm($i) 1.0e-4

        # Output material information to data file
        puts $materialProperties "
        Material ID: $i
        rho: $rho($i) 
        G: $G($i)
        K: $bulk($i)
        Phi: $phi($i)
        GammaPeak: $gammaPeak
        refPress: $refPress($i)
        pressCoeff: $pressCoeff($i)
        phaseAng: $phaseAng($i)
        contract1: $contract1($i)
        contract3: $contract3($i)
        dilate1: $dilate1($i)
        dilate3: $dilate3($i)
        "

} else {
    
	# PIMY Model
    nDMaterial PressureIndependMultiYield $i 2 $rho($i) $G($i) $bulk($i) $cohesion($i) $gammaPeak 0.0 \
    $refPress($i) $pressCoeff($i) 20

    set thick($i) 1.0
    set xWgt($i)  [expr $g*sin($slope)]
    set yWgt($i)  [expr $g*cos($slope)]
    set porosity($i) [expr $voidR($i) / (1 + $voidR($i))]
    set uBulk($i) [expr $bulkWater/$porosity($i)]
    set hPerm($i) 1.0e-4
    set vPerm($i) 1.0e-4

    # Output material information to data file
    puts $materialProperties "
    Material ID: $i
    rho: $rho($i) 
    G: $G($i)
    K: $bulk($i)
    coheison: $cohesion($i)
    Phi: 0.0
    GammaPeak: $gammaPeak
    refPress: $refPress($i)
    pressCoeff: $pressCoeff($i) 
    "

  }
}

close $materialProperties

# Define linear elastic material for "bedrock"
nDMaterial ElasticIsotropic 0 $E(0) 0.3 $rho(0)

set thick(0) 1.0
set xWgt(0)  [expr $g*sin($slope)]
set yWgt(0)  [expr $g*cos($slope)]
set porosity(0) [expr $voidR(0) / (1 + $voidR(0))]
set uBulk(0) [expr $bulkWater/$porosity(0)]
set hPerm(0) 1.0e-4
set vPerm(0) 1.0e-4

puts "Finished creating all soil materials..."

#-----------------------------------------------------------------------------------------
#  6. CREATE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------

# Define the top element number of each layer
set elemCount 0.0
set layerTopEleNum(-1) 0
for {set i 0} {$i <= $numLayers} {incr i 1} {
    set layerTopEleNum($i) [expr $elemCount + $nElemY($i)]
    set elemCount [expr $elemCount + $nElemY($i)]
}

for {set j 1} {$j <= $nElemT} {incr j 1} {

    set nI  [expr 6*$j - 5]
    set nJ  [expr $nI + 2]
    set nK  [expr $nI + 8]
    set nL  [expr $nI + 6]
    set nM  [expr $nI + 1]
    set nN  [expr $nI + 5]
    set nP  [expr $nI + 7]
    set nQ  [expr $nI + 3]
    set nR  [expr $nI + 4]

    set lowerBound 0.0
    for {set i 0} {$i <= $numLayers} {incr i 1} {

        if {$j <= $layerTopEleNum($i) && $j > $layerTopEleNum([expr $i - 1])} {

            # permeabilities are initially set at 10.0 m/s for gravity analysis, values are updated post-gravity
            element 9_4_QuadUP $j $nI $nJ $nK $nL $nM $nN $nP $nQ $nR \
            $thick($i) $i $uBulk($i) 1.0 10.0 10.0 $xWgt($i) $yWgt($i)

            puts "element 9_4_QuadUP $j $nI $nJ $nK $nL $nM $nN $nP $nQ $nR $thick($i) $i $uBulk($i) 1.0 1.0 1.0 $xWgt($i) $yWgt($i)"

    }
    set lowerBound $layerBound($i)
}
}
puts "Finished creating all soil elements..."

# print all node info to text file (if statement so it only does it once when batching)
if {[file exists "$saveDir/$siteModelID.eleInfo.txt"] == "0"} {
     print $saveDir/$siteModelID.eleInfo.txt -ele
}
#-----------------------------------------------------------------------------------------
#  7. LYSMER DASHPOT
#-----------------------------------------------------------------------------------------

# Define dashpot nodes
set dashF [expr $nNodeT+1]
set dashS [expr $nNodeT+2]
#    nodeTag x   y
node $dashF  0.0 0.0
node $dashS  0.0 0.0
# Define fixities for dashpot nodes
#   nodeTag dof1 dof2
fix $dashF  1    1
fix $dashS  0    1

# Define equal DOF for dashpot and base soil node
#        rNodeTag cNodeTag dof1
equalDOF 1        $dashS   1
puts "Finished creating dashpot nodes and boundary conditions..."

# Define dashpot material
set colArea       [expr $sElemX*$thick(0)]
set dashpotCoeff  [expr $rockVS*$rockDen]
#                         matTag              C: damping coefficient       alpha (1=linear damping)
uniaxialMaterial Viscous [expr $numLayers+1] [expr $dashpotCoeff*$colArea] 1

# Define dashpot element
#                   eleTag           iNode  jNode        matTag                   direction (1=X)
element zeroLength [expr $nElemT+1]  $dashF $dashS -mat [expr $numLayers+1]  -dir 1
puts "Finished creating dashpot material and element..."

#-----------------------------------------------------------------------------------------
#  8. CREATE GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------

# Create list for pore pressure nodes
set nodeList3 {}
set channel [open "$saveDir/$siteModelID.ppNodesInfo.dat" r]
set count 0;
foreach line [split [read -nonewline $channel] \n] {
    set count [expr $count+1];
    set lineData($count) $line
    set nodeNumber [lindex $lineData($count) 0]
    lappend nodeList3 $nodeNumber
}
set nodeList3 [lsort -integer $nodeList3]
puts "pore pressure nodes are: $nodeList3"
close $channel

# record nodal displacment, acceleration, and porepressure
#eval "recorder Node -file Gdisplacement.out -time -node $nodeList3 -dof 1 2  disp"
#eval "recorder Node -file Gacceleration.out -time -node $nodeList3 -dof 1 2  accel"
#eval "recorder Node -file GporePressure.out -time -node $nodeList3 -dof 3 vel"
# record elemental stress and strain (files are names to reflect GiD gp numbering)
#recorder Element -file Gstress1.out   -time  -eleRange 1 $nElemT  material 1 stress
#recorder Element -file Gstress2.out   -time  -eleRange 1 $nElemT  material 2 stress
#recorder Element -file Gstress3.out   -time  -eleRange 1 $nElemT  material 3 stress
#recorder Element -file Gstress4.out   -time  -eleRange 1 $nElemT  material 4 stress
#recorder Element -file Gstress9.out   -time  -eleRange 1 $nElemT  material 9 stress
#recorder Element -file Gstrain1.out   -time  -eleRange 1 $nElemT  material 1 strain
#recorder Element -file Gstrain2.out   -time  -eleRange 1 $nElemT  material 2 strain
#recorder Element -file Gstrain3.out   -time  -eleRange 1 $nElemT  material 3 strain
#recorder Element -file Gstrain4.out   -time  -eleRange 1 $nElemT  material 4 strain
#recorder Element -file Gstrain9.out   -time  -eleRange 1 $nElemT  material 9 strain
puts "Finished creating gravity recorders..."

#-----------------------------------------------------------------------------------------
#  9. DEFINE ANALYSIS PARAMETERS
#-----------------------------------------------------------------------------------------

#---GROUND MOTION PARAMETERS

# Number of steps in ground motion record. 
# Retrieved from GM file name (exterior input)
set nSteps $npts
puts "duration = [expr $nSteps*$dt] sec"

#---RAYLEIGH DAMPING PARAMETERS
set pi      3.141592654
# Damping ratio
set damp    0.05
# Lower frequency
#set omega1  [expr 2*$pi*Tsite]
set omega1  [expr 2*$pi*0.2]
# Upper frequency
#set omega2  [expr 2*$pi*5*Tsite]
set omega2  [expr 2*$pi*20]
# Damping coefficients
set a0      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set a1      [expr 2*$damp/($omega1 + $omega2)]
puts "damping coefficients: a_0 = $a0;  a_1 = $a1"

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
#STABILITY CHECK NOT REQUIRED FOR IMPLICIT (NEWMARK) ANALYSIS
# maximum shear wave velocity (m/s)
#set vsMax $Vs(0)
#for {set i 1} {$i <= $numLayers} {incr i 1} {
#    if {$Vs($i) > $vsMax} {
#        set vsMax $Vs($i)
#}
#}
# duration of ground motion (s)
#set duration    [expr $motionDT*$motionSteps]
# minimum element size
#set minSize $sElemY(0)
#for {set i 1} {$i <= $numLayers} {incr i 1} {
#    if {$sElemY($i) < $minSize} {
#        set minSize $sElemY($i)
#}
#}
#
# trial analysis time step
#set kTrial      [expr $minSize/(pow($vsMax,0.5))]
# define time step and number of steps for analysis
#if { $motionDT <= $kTrial } {
#    set nSteps  $motionSteps
#    set dt      $motionDT
#} else {
#    set nSteps  [expr int(floor($duration/$kTrial)+1)]
#    set dt      [expr $duration/$nSteps]
#}
puts "number of steps in analysis: $nSteps"
puts "analysis time step: $dt"

#---ANALYSIS PARAMETERS

# Newmark parameters
# Average acceleration method:
# Unconditionally stable, without numerical damping
set gamma           0.50
set beta            0.25

#-----------------------------------------------------------------------------------------
#  10. GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

# Update materials to ensure elastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 0
}

# Create analysis
constraints Penalty 1.e14 1.e14
test        NormDispIncr 1e-6 35 1
algorithm   KrylovNewton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta
analysis    Transient

set startT  [clock seconds]

# Elastic gravity loading
# Large time steps are used to damp out the waves generated by loading
#       numincr dt
analyze 10      5.0e2
puts "Finished with elastic gravity analysis..."

# Update materials to consider plastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 1
}

# Plastic gravity loading
#       numincr dt
analyze 40      5.0e-2
puts "Finished with plastic gravity analysis..."

#-----------------------------------------------------------------------------------------
#  11. UPDATE ELEMENT PERMEABILITY VALUES FOR POST-GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

# If excess pore pressure generation is not allowed (i.e., allowPWP = No),
# permeabilites are left high for dynamic analysis

if {$allowPWP == Yes} {

    # Choose base number for parameter IDs which is higer than other tags used in analysis
    set ctr 10000.0
    # Loop over elements to define parameter IDs 
    for {set i 1} {$i<=$nElemT} {incr i 1} {
        parameter [expr int($ctr+1.0)] element $i vPerm
        parameter [expr int($ctr+2.0)] element $i hPerm
        set ctr [expr $ctr+2.0]
}

# Update permeability parameters for each element using parameter IDs
set ctr 10000.0
for {set j 1} {$j <= $nElemT} {incr j 1} {

    set lowerBound 0.0
    for {set i 0} {$i <= $numLayers} {incr i 1} {

        if {$j <= $layerTopEleNum($i) && $j > $layerTopEleNum([expr $i - 1])} {
            updateParameter [expr int($ctr+1.0)] $vPerm($i)
            updateParameter [expr int($ctr+2.0)] $hPerm($i)
            puts "$j  updateParameter [expr int($ctr+1.0)] $vPerm($i)"

        }
        set lowerBound $layerBound($i)
    }
    set ctr [expr $ctr+2.0]
}
puts "Finished updating permeabilities for dynamic analysis..."
}

#-----------------------------------------------------------------------------------------
#  12. CREATE POST-GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------

# Create list for pore pressure nodes
set leftPPNodeList {}
set channel1 [open "$saveDir/$siteModelID.LeftPPNodesInfo.dat" r]
set count1 0;
foreach line [split [read -nonewline $channel1] \n] {
    set count1 [expr $count1+1];
    set lineData($count1) $line
    set nodeNumber [lindex $lineData($count1) 0]
    lappend leftPPNodeList $nodeNumber
    puts "node: $nodeNumber"
}
set leftPPNodeList [lsort -integer $leftPPNodeList]
puts "leftPPNodes are: $leftPPNodeList"
close $channel1

# Reset time and analysis
setTime 0.0
wipeAnalysis
remove recorders

# Recorder time step
set recDT  [expr 10*$dt]

# Record nodal displacment, acceleration, and porepressure
#eval "recorder Node -file $saveDir/$analID.disp.out -time -dt $dt -node $leftPPNodeList -dof 1 2  disp"
#eval "recorder Node -file $saveDir/$analID.vel.out -time -dt $dt -node $leftPPNodeList -dof 1 2  vel"
#eval "recorder Node -file $saveDir/$analID${dash}acc.out -time -node $leftPPNodeList -dof 1 2 accel"
#recorder Node -file $saveDir/$analID${dash}accHorSur.out -time -node $nNodeT -dof 1 accel
# record horizontal acceleration at the top node
eval "recorder Node -file $saveDir/$analID${dash}accHorSur.out -time -node $nNodeT -dof 1 accel"

# record elemental stress and strain
#recorder Element -file $saveDir/$analID.stress1.out   -time  -eleRange 1 $nElemT  material 1 stress 3
#recorder Element -file $saveDir/$analID.stress2.out   -time  -eleRange 1 $nElemT  material 2 stress 3
#recorder Element -file $saveDir/$analID.stress3.out   -time  -eleRange 1 $nElemT  material 3 stress 3
#recorder Element -file $saveDir/$analID.stress4.out   -time  -eleRange 1 $nElemT  material 4 stress 3
#recorder Element -file $saveDir/$analID.stress5.out   -time  -eleRange 1 $nElemT  material 5 stress 3
#recorder Element -file $saveDir/$analID.stress6.out   -time  -eleRange 1 $nElemT  material 6 stress 3
#recorder Element -file $saveDir/$analID.stress7.out   -time  -eleRange 1 $nElemT  material 7 stress 3
#recorder Element -file $saveDir/$analID.stress8.out   -time  -eleRange 1 $nElemT  material 8 stress 3
#recorder Element -file $saveDir/$analID${dash}stress9.out   -time  -eleRange 1 $nElemT  material 9 stress 3
#recorder Element -file $saveDir/$analID.strain1.out   -time  -eleRange 1 $nElemT  material 1 strain
#recorder Element -file $saveDir/$analID.strain2.out   -time  -eleRange 1 $nElemT  material 2 strain
#recorder Element -file $saveDir/$analID.strain3.out   -time  -eleRange 1 $nElemT  material 3 strain
#recorder Element -file $saveDir/$analID.strain4.out   -time  -eleRange 1 $nElemT  material 4 strain
#recorder Element -file $saveDir/$analID.strain5.out   -time  -eleRange 1 $nElemT  material 5 strain
#recorder Element -file $saveDir/$analID.strain6.out   -time  -eleRange 1 $nElemT  material 6 strain
#recorder Element -file $saveDir/$analID.strain7.out   -time  -eleRange 1 $nElemT  material 7 strain
#recorder Element -file $saveDir/$analID.strain8.out   -time  -eleRange 1 $nElemT  material 8 strain
#recorder Element -file $saveDir/$analID${dash}strain9.out   -time  -eleRange 1 $nElemT  material 9 strain
puts "Finished creating all recorders..."


#-----------------------------------------------------------------------------------------
#  13. DYNAMIC ANALYSIS
#-----------------------------------------------------------------------------------------

# Create a 2D model with 3 degrees of freedom
model BasicBuilder -ndm 2 -ndf 3

# Define constant scaling factor for applied velocity
set cFactor [expr $colArea*$dashpotCoeff]

# Define velocity time history file
set velocityFile $gMotionPath

# Timeseries object for force history
set mSeries "Path -dt $dt -filePath $velocityFile -factor $cFactor"

# Loading object
#             patternTag timeSeries
pattern Plain 13         $mSeries {
#        nodeTag dof1 dof2 dof3
    load 1       1.0  0.0  0.0
}
puts "Dynamic loading created..."

# Create analysis
constraints Penalty 1.e16 1.e16
test        NormDispIncr 1.0e-5 35 1
algorithm   KrylovNewton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta
rayleigh    $a0 $a1 0.0 0.0
analysis    Transient

# Perform analysis with timestep reduction loop
#               numincr	dt
set ok [analyze $nSteps $dt]

# If analysis fails, reduce timestep and continue with analysis
if {$ok != 0} {
    puts "did not converge, reducing time step"
    set curTime  [getTime]
    set mTime $curTime
    puts "curTime: $curTime"
    set curStep  [expr $curTime/$dt]
    puts "curStep: $curStep"
    set rStep  [expr ($nSteps-$curStep)*2.0]
    set remStep  [expr int(($nSteps-$curStep)*2.0)]
    puts "remStep: $remStep"
    set dt       [expr $dt/2.0]
    puts "dt: $dt"
    puts "Ground Motion is: $gMotionName"

    set ok [analyze  $remStep  $dt]

    # If analysis fails again, reduce timestep and continue with analysis
    if {$ok != 0} {
        puts "did not converge, reducing time step"
        set curTime  [getTime]
        puts "curTime: $curTime"
        set curStep  [expr ($curTime-$mTime)/$dt]
        puts "curStep: $curStep"
        set remStep  [expr int(($rStep-$curStep)*2.0)]
        puts "remStep: $remStep"
        set dt       [expr $dt/2.0]
        puts "dt: $dt"
        puts "Ground Motion is: $gMotionName"

        analyze  $remStep  $dt
}
}
set endT    [clock seconds]
puts "Finished with dynamic analysis..."
puts "Analysis execution time: [expr $endT-$startT] seconds"

wipe
