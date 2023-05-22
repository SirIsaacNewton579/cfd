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

### 1.
(pL,rhoL,uL,pR,rhoR,uR)=(1,1,-0.,0.1,0.125,-0.);  

The gif of the solution:  

rho:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/rhoo.gif )     

u:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/uo.gif )     

p:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/po.gif )    

### 2.
(pL,rhoL,uL,pR,rhoR,uR)=(0.1,0.125,0,1,1,0);    

left-going shock wave, right-going expansion wave;  

The gif of the solution:  

rho:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/rholp.gif )     

u:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/ulp.gif )     

p:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/plp.gif )    

### 3.
(pL,rhoL,uL,pR,rhoR,uR)=(1,1,-1,0.5,0.5,1);  

left-going shock expansion wave, right-going expansion wave;  

The gif of the solution:  

rho:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/rhopp.gif )     

u:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/upp.gif )     

p:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/ppp.gif )   

### 4.
(pL,rhoL,uL,pR,rhoR,uR)=(1,1,-10,1,1,10);  

left-going shock expansion wave, right-going expansion wave, vacuum in the middle;  

The gif of the solution:  

rho:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/rhoppv.gif )     

u:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/uppv.gif )     

p:  

![image](https://github.com/SirIsaacNewton579/cfd/blob/master/pppv.gif )   
