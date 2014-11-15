function [a, b, c, d, e, f, g, h] = ChemicalEquationBalance(oxMassFlow)

%Set up initial conditions
a = 10; %oxMassFlow; 
b = 1/((0.6/0.4)*(26.98/352.52)); %# of mols for aluminum on reactant side
h = a;

reactants = [2 1 3 0 0; 1 0 0 25 0; 0 2 0 52 0; 0 0 2 0 1; 0 0 0 b -1];
products = [a; 25; 52; b; 0];

x = reactants\products;

c = x(1);
d = x(2);
e = x(3);
f = x(4);
g = x(5);

