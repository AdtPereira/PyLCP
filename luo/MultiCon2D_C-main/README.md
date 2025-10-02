# 2D capacitance calculation method for round conductors with insulation layer

This method is based on analytical solution under static electric  field for round condcutors.

This code is built in MATLAB 2019b

**Example**

There are 2 conductors with 1mm radius and 1mm thickness insulation layer, permittivity is 2.2. They locate at (0,0) and (4.4e-3,0), unit m, respectively.

To calculate the charges on the surface of conductors with opposite 0.5 V, and the order is set to 6, commend is listed as follow

x = [0;4.4e-3]; y = [0;0]; a = 1e-3*ones(2,1); b = 2e-3*ones(2,1); V = [0.5,-0.5];

Q = MultiConMatrix(x,y,a,'b',b,'V',V,'Nord',6,'kind','C2','er2',2.2);

Result Q is 

2.9256e-11

-2.9256e-11

Result from COMSOL is

2.9188E-11

-2.9188E-11

**Author:**

Tianming Luo, TU Delft, HV group

https://www.tudelft.nl/ewi/over-de-faculteit/afdelingen/electrical-sustainable-energy/high-voltage-technologies/research/design-method-for-high-power-medium-frequency-transformers