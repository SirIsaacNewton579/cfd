# cfd
Exact Solution of Riemann problems;
LeftState=(rhoL,uL,pL);
RightState=(rhoR,uR,pR);

## how to use

To solve different situations, Just make some changes on line 236 of "cfdex1.cpp": 
SetState(pL,rhoL,uL,pR,rhoR,uR); 
where pL,rhoL,uL,pR,rhoR,uR are the initial state;

## examples

Here is the example.
(pL,rhoL,uL,pR,rhoR,uR)=(1,1,-0.,0.1,0.125,-0.);
The gif of the solution:
rho:
![image](https://github.com/SirIsaacNewton579/cfd/blob/master/rho.gif )   
u:
![image](https://github.com/SirIsaacNewton579/cfd/blob/master/u.gif )   
p:
![image](https://github.com/SirIsaacNewton579/cfd/blob/master/p.gif )   
