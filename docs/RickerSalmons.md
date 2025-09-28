---
layout: page
author: Jordi Villà-Freixa
title: RickerSalmons
permalink: /RickerSalmons
---
<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# **1.1.1. Models Unidimensionals. Model de Ricker**
## **Cas: el salmó roig al nord del pacífic**

![image_0.png](RickerSalmons_media/image_0.png)


El salmó roig, als rius del nord del Pacífic, fa una posta cada 4 anys. Si en tres períodes consecutius, $t=0,4,8$ , s'observen poblacions de 0.325, 0.431 i 0.529 milions de salmons, respectivament, és possible trobar els valors de $R$ i $K$ que permetin calcular l'evolució de la població de salmons per un model de Ricker?


**Resolució amb Matlab del sistema d'equacions de l'exemple del Model de Ricker dels salmons als rius del nord del Pacífic**


*Es tracta d'un sistema de dues equacions amb dues incògnites: r i K ( i tot i que no són equacions lineals, el Matlab no té problemes en trobar una solució)*

```matlab
% entrem la instrucció per a resoldre sistemes d'equacions que figura al
% document "Introducció al Matlab"

syms r K, solucio=solve([0.431==0.325*exp(r*(1-0.325/K)),0.529==0.431*exp(r*(1-0.431/K))],[r,K])
```

```matlabTextOutput
Warning: Unable to solve symbolically. Returning a numeric solution using vpasolve.
solucio = struct with fields:
    r: 0.51960209604121860724470879574985
    K: 0.71157617909621780281959432221705

```

Tot i que el programa ens avisa que no pot resoldre el sistema analíticament, i ens proposa de resoldre'l numèricament a partir d'una instrucció diferent del "solve", que és "vpasolve", fixeu\-vos que ja ens ha donat la solució que buscàvem. 


Però si no n'esteu segurs i voleu repetir\-ho fent servir aquesta nova instrucció "vpasolve", ho fem a continuació:


Això sí, reanomeno les incògnites r i K com r1 i K1, per tal que el programa no es confongui. També reanomeno la solució com a solucio1

```matlab
syms r1 K1, solucio1=vpasolve([0.431==0.325*exp(r1*(1-0.325/K1)),0.529==0.431*exp(r1*(1-0.431/K1))],[r1,K1])
```

```matlabTextOutput
solucio1 = struct with fields:
    r1: 0.51960209604121860724470879574985
    K1: 0.71157617909621780281959432221705

```

En aquest cas de solució numèrica, no ens cal ni demanar\-li que ens mostri la solució com hauríem de fer en una solució analítica on caldria la instrucció següent:

```matlab
[solucio1.r1,solucio1.K1]
```
ans = 
 $\displaystyle \left(\begin{array}{cc} 0.51960209604121860724470879574985 & 0.71157617909621780281959432221705 \end{array}\right)$
 

I si ho preferim, aquests resultats de les variables r1 i K1 només amb 4 decimals, la instrucció "single" ens ho permet:

```matlab
single(ans)
```

```matlabTextOutput
ans = 1x2 single row vector
    0.5196    0.7116

```
