---
layout: page
author: Jordi Villà-Freixa
title: DiscretParasitHoste
permalink: /DiscretParasitHoste
---
<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Model discret: Dinàmica paràsit\-hoste
\matlabtableofcontents

Considerem un model discret on les poblacions del **paràsit** $P_t$ i l'**hoste** $H_t$ evolucionen al llarg del temps en intervals $t$ . Considerem dos models.

## Model de Nicholson\-Bailey
 $$ \begin{array}{l} H_{k+1} =RH_k e^{-aP_k } \newline P_{k+1} =SH_k \left(1-e^{-aP_k } \right) \end{array} $$ 

on 

-  $R$ : constant de creixement dels hostes; 
-  $a$ : eficiència dels paràsits en cercar hostes; i 
-  $S$ : promig d'ous viables dels paràsits per a cada hoste infectat. 

## Model bionomial negatiu (Griffiths\-May)

Usant la mateixa notació que en el cas anterior, aquest model representa les poblacions d'hostes i paràsits seguint:

 $$ \begin{array}{l} H_{k+1} =RH_k {\left(1+\frac{aP_k }{m}\right)}^{-m} \newline P_{k+1} =SH_k \left(1-{\left(1+\frac{aP_k }{m}\right)}^{-m} \right) \end{array} $$ 

on apareix el nou paràmetre $m$ per controlar l'efecte de l'eficiència del paràsit.

## Qüestions
1.  Construeix gràfiques per al model de Nicholson\-Bailey usant $R=1.5$ , $S=1$ i $a=0.023$ per als casos: a) $H_0 =30$ i $P_0 =20$ i b) $H_0 =20$ i $P_0 =20$ (Figura 2.8 dels apunts de models discrets). Com interpretes la diferència entre les dues?
2. Usant els mateixos paràmetres anteriors i també $m=0.5$ , repeteix els dos gràfics anteriors per al model binomial negatiu. Quina conclusió n'extreus?
3. Construeix un gràfic de fase (població dels paràsits en funció de la població dels hostes) per als diferents casos. Què observes?
4. Comprova gràficament l'estabilitat dels dos punts d'equilibri trobats a l'exemple 2.55 dels apunts de models discrets.

## Respostes

**Algunes instruccions preliminars:**

```matlab
startup
```

Codi MATLAB per simular el model discret de l'exemple 2.56 dels apunts

```matlab
% Paràmetres del model
tmax = 10;  % Temps màxim de simulació
x0 = 20/3*log(2)+0.01;  % Població inicial d'hostes
y0 = 10*log(2);  % Població inicial de paràsits

% Inicialització dels vectors de població
x = zeros(1, tmax);
y = zeros(1, tmax);
x(1) = x0;
y(1) = y0;

% Simulació del sistema discret
for t = 1:tmax-1
    x(t+1) = 2*x(t)*exp(-0.1*y(t));
    y(t+1) = 3*x(t)*(1-exp(-0.1*y(t)));
end

% Representació gràfica de les poblacions
figure;
plot(1:tmax, x, 'b', 1:tmax, y, 'r');
title('Dinàmica de paràsit-hoste en temps discret');
xlabel('Temps');
ylabel('Població');
legend('Hostes', 'Paràsits');

% Gràfic en fase (Població d'hostes vs paràsits)
figure;
plot(x, y, 'k', 'LineWidth', 1.5);
title('Gràfic en fase: Població d\''Hostes vs Paràsits');
xlabel('Població d\''Hostes');
ylabel('Població de Paràsits');
grid on;
```

Codi MATLAB per simular el model discret de l'exemple 2.55 dels apunts

```matlab
% Paràmetres del model
tmax = 100;  % Temps màxim de simulació
x0 = 0.1;  % Població inicial d'hostes
y0 = .01;  % Població inicial de paràsits

% Inicialització dels vectors de població
x = zeros(1, tmax);
y = zeros(1, tmax);
x(1) = x0;
y(1) = y0;

% Simulació del sistema discret
for t = 1:tmax-1
    x(t+1) = y(t);
    y(t+1) = 0.5*x(t)+y(t)-(y(t)).^2;
end

% Representació gràfica de les poblacions
figure;
plot(1:tmax, x, 'b', 1:tmax, y, 'r');
title('Dinàmica de paràsit-hoste en temps discret');
xlabel('Temps');
ylabel('Població');
legend('Hostes', 'Paràsits');

% Gràfic en fase (Població d'hostes vs paràsits)
figure;
plot(x, y, 'k', 'LineWidth', 1.5);
title('Gràfic en fase: Població d\''Hostes vs Paràsits');
xlabel('Població d\''Hostes');
ylabel('Població de Paràsits');
grid on;
```
