# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 13:03:27 2021

@author: Sarah Chancellor, UCCS Undergraduate Mechanical Engineering Department. 
AIAA propulsion proposal and calculations


All units are in kg for mass, seconds for time, pressure is in Pa, Kelvin in Temperature, J or enegery, and m for distance 
SF is the scaling factor with respect to the SpaceX Raptor rocket engine 
""" 

"""
max G's considering that space travel is taxing on the human body, the force 
experienced on take-off should be around 1-2 earth G's. Astronauts experience 
around 3 exiting earth atmosphere, not as much force is needed on Mars. 
Constant acceleration will be established after the inital first moments.
Calculations are based on constant acceleration.  
"""
import math 
import CoolProp

maxalt= 450* 1000
fuel= 1000
g= 3.711
m= 2000 
minthrust1= g*m  
print(minthrust1)
SF= 0.10
M=1  
"""
Nozzle calculations: 
    i subscript stands for inlet or throat of nozzle, e stands for the exit
    c is for the chamber which comes before the throat. m for methane and o for Oxygen
    CH4 + 2O2 = CO2 + 2H2O
  
"""
Mmo2= 32
Mmch4= (4*1.008) +12.01
Mh20= (18.015)
Mco2 = 12.01 + 32
minnozzleA= pow(2,2)* math.pi
mnd= (pow(SF,0.483))*2.98704
maxnozzleA= pow(mnd,2) * math.pi 
Tc= 749
Pc= pow(SF,-0.4968)*300*100000
Pi= Pc
cpm= CoolProp.CoolProp.PropsSI('C', 'P', Pi, 'T',Tc, 'Methane')
cpo= CoolProp.CoolProp.PropsSI('C', 'P', Pi, 'T',Tc, 'Oxygen')
cvm= CoolProp.CoolProp.PropsSI('O', 'P', Pi, 'T',Tc, 'Methane')
cvo= CoolProp.CoolProp.PropsSI('O', 'P', Pi, 'T',Tc, 'Oxygen')
k1r= ((Mmch4/(Mmch4+Mmo2))*cvm) + ((Mmo2/(Mmch4+Mmo2))*cvo)
k2r= ((Mmch4/(Mmch4+Mmo2))*cpm) + ((Mmo2/(Mmch4+Mmo2))*cpo)
kr= k2r/k1r
Te= Tc* pow((1+ (0.5*(kr-1)*pow(M,2))), (-1))
cpw= CoolProp.CoolProp.PropsSI('C', 'P', Pi, 'T',Te, 'Water')
cpcd= CoolProp.CoolProp.PropsSI('C', 'P', Pi, 'T',Te, 'CarbonDioxide')
cvw= CoolProp.CoolProp.PropsSI('O', 'P', Pi, 'T',Te, 'Water')
cvcd= CoolProp.CoolProp.PropsSI('O', 'P', Pi, 'T',Te, 'CarbonDioxide')
k1p= ((Mh20/(Mh20+ Mh20+ Mco2))*cvw) + ((Mco2/(Mh20+Mh20+Mco2))*cvcd)
k2p= ((Mh20/(Mh20+Mh20+Mco2))*cpw) + ((Mco2/(Mh20+Mh20+Mco2))*cpcd)
kp= k2p/k1p
vi= 0
rw= CoolProp.CoolProp.PropsSI('gas_constant', 'P', Pi, 'T',Tc, 'Water')
rcd= CoolProp.CoolProp.PropsSI('gas_constant', 'P', Pi, 'T',Tc, 'CarbonDioxide')
R= (((2*Mh20)/((2*Mh20)+Mco2))*(rw)) + ((Mco2/((2*Mh20)+Mmo2))*(rcd))
"In comparison to the exit velocity, the inlet velocity can be approximated as close to zero. "
Pe= Pc * pow((1+ (0.5*(kp-1)*pow(M,2))), (-kp/(kp-1)))
ve2= (2*Pc*Tc*R*kp)/(kp-1)
ve1= 1- pow((Pe/Pi),((kp-1)/kp)) 
        
vexit= pow((ve1*ve2),0.5)
  
romi =CoolProp.CoolProp.PropsSI('D','P', Pi, 'T', Tc, 'Methane')
rooi = CoolProp.CoolProp.PropsSI('D', 'P', Pi, 'T',Tc, 'Oxygen')
ro= (Mmch4/(Mmch4+Mmo2))*romi + (Mmo2/(Mmch4+Mmo2))*rooi 
mdotfuel= (1/1000) * (maxnozzleA*Pc*(pow(Tc, -0.5))) *(pow((kp/R),0.5)) *M * pow((1+ (0.5*(kp-1))*pow(M,2)), ((-kp+1)/(2*(kp-1))))
Thrust= ((mdotfuel)*(vexit)) + ((Pe-Pi)*maxnozzleA)
t= 60*2
fuel= t*mdotfuel
minthrust= (m+fuel)*g
 
print('minthrust= '+ str(minthrust) + ' KN')
print( 'thrust= '+str(Thrust) + ' KN')
print('mass flow rate= ' +str(mdotfuel)+ ' kg/s')
print('Max total fuel mass= ' + str(fuel) + ' kg. Note: this number is if the rocket is going full power the entire time wich will not be the case. Safety measure is inferred. the necessary thrust will decrease as the mass increases and the MAV goes into space.')
