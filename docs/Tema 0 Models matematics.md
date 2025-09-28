---
layout: page
author: Jordi Villà-Freixa
title: "Tema 0. Models Matemàtics en Biologia"
permalink: /Tema0
---

::: frame
:::

::: frame
### Índex
:::

# Prèvia

::::: frame
### Models matemàtics en biologia

:::: columns
::: center
![image](mathmodels.jpg){height="70%"}
:::

El procés de modelització en les ciències de la vida. Les activitats
principals implicades en aquest procediment són l'observació seguida de
la modelització matemàtica; simulació, anàlisi, optimització i retorn a
l'observació. En aquest cicle, el model matemàtic ocupa, just després
del sistema real, la posició central. Extret de
[@torres_mathematical_2015].
::::
:::::

::: frame
Material

### Referències

El material d'aquestes presentacions està basat en anteriors
presentacions i apunts d'altres professors [@corbera2019] de la
UVic-UCC, pàgines web diverses (normalment enllaçades des del text),
així com monografies [@otto_biologists_2007].
:::

::: frame
Introducció al curs: Matemàtiques aplicades a la Biologia

-   Objectiu principal: millorar la capacitat d'entendre i construir
    models matemàtics en biologia.

-   Les matemàtiques són com una llengua estrangera per a molts biòlegs:
    cal temps i pràctica per dominar-les.

-   L'important no és memoritzar fórmules, sinó entendre-les i saber on
    buscar quan no coneixes la solució.

-   Matemàtiques utilitzades en biologia: àlgebra, càlcul bàsic, àlgebra
    de matrius, teoria de probabilitats, i estadística.

-   Avantatge de l'era dels ordinadors: no cal ser un expert en
    matemàtiques per resoldre models complexos.

-   Ús del programari Matlab per realitzar càlculs matemàtics avançats,
    permetent centrar-se en la formulació i anàlisi de models.
:::

# Posem un exemple

::: frame
Exemple: Població de Porcs Senglar Per exemple, si volem estudiar la
població futura de porcs senglars en un hàbitat, podem utilitzar un
model matemàtic en lloc de realitzar experiments físics, cosa que és més
barata i eficient.\
Així doncs, plantejarem models matemàtics i utilitzarem eines per a
resoldre'ls.
:::

# Models matemàtics

:::: frame
Definició de Model Matemàtic Un model matemàtic és una representació
d'un procés real mitjançant equacions. Per construir un model matemàtic
s'ha de:

-   Abstreure la realitat per identificar els elements importants.

-   Definir variables significatives i interpretar les seves relacions.

-   Simplificar i comprovar la validesa del model.

::: epigraph
"Everything should be made as simple as possible, but not simpler\"

Albert Einstein
:::
::::

::: frame
Mètode Científic en la construcció de models Els mètodes matemàtics es
basen en el mètode científic:

1.  Anàlisi inicial de la situació: cal esbrinar quines són les
    qüestions més rellevants, les quals ajudin a simplificar el problema
    que s'estudia

2.  Interpretació del funcionament qualitatiu del sistema: es dissenya
    un model molt simple que descrigui el comportament del sistema

3.  Descripció quantitativa del model: cal definir les variables
    significatives del sistema i les seves interrelacions (hipòtesis del
    model)

4.  Escriptura de les equacions que regeixen el comportament del sistema

5.  Comprovació de la validesa del model per a algun cas conegut

6.  Reanàlisi de les hipòtesis inicialment plantejades si la comprovació
    anterior no és bona (tornant al pas número 3) o, en cas contrari,
    acceptació del model dissenyat

7.  Utilització del model per a fer prediccions del comportament del
    sistema
:::

::: frame
Classificació de Models Matemàtics Els models matemàtics es classifiquen
segons l'evolució de les variables en el **temps**:

-   Models discrets: Equacions recursives. Sovint matricials.

-   Models continus: Equacions diferencials (ordinàries -ODE- o parcials
    -PDE-).

Els models discrets s'estudien a la primera part del curs i els models
continus a la segona.
:::

:::: frame
::: center
[![image](Discrete-vs-continuum.jpg){width="70%"}](https://www.zoology.ubc.ca/~bio301/Bio301/Lectures/Lecture7/Overheads.html)
:::
::::

::: frame
Set passos per modelar un problema biològic[@otto_biologists_2007] **Pas
1: Formular la pregunta**

-   Què vols saber?

-   Descriu el model en forma de pregunta.

-   Simplifica la pregunta.

-   Comença amb la descripció més simple i raonable del problema.

**Pas 2: Determinar els ingredients bàsics**

-   Defineix les variables del model i les seves restriccions.

-   Descriu les interaccions entre variables.

-   Decideix si el temps serà discret o continu, i estableix l'escala
    temporal.

-   Defineix els paràmetres del model i les seves restriccions.

**Pas 3: Descriure qualitativament el sistema biològic**

-   Dibuixa un diagrama de cicle de vida o de flux per descriure els
    canvis en les variables.

-   Si hi ha molts esdeveniments possibles, construeix una taula amb els
    resultats de cada esdeveniment.
:::

::: frame
Set passos per modelar un problema biològic (cont.) **Pas 4: Descriure
quantitativament el sistema biològic**

-   Escriu les equacions basant-te en els diagrames i taules.

-   Fes comprovacions: es compleixen les restriccions a mesura que passa
    el temps? Els resultats del model poden abordar la pregunta?

**Pas 5: Analitzar les equacions**

-   Simula i dibuixa gràfics dels canvis en el sistema.

-   Realitza les anàlisis apropiades i assegura't que aborden el
    problema.

**Pas 6: Verificacions**

-   Comprova els resultats amb dades o casos especials coneguts.

-   Considera alternatives al model més simple i repeteix els passos
    anteriors.

**Pas 7: Relacionar els resultats amb la pregunta**

-   Els resultats responen a la pregunta biològica?

-   Interpreta els resultats verbalment i descriu nous coneixements
    sobre el procés biològic.

-   Descriu possibles experiments.
:::

# Apliquem els models a un cas concret

::: frame
Exemple 1a: Model Discret de l'evolució d'un cultiu ceular

Pas 1:

:   El nombre de cèules creix; a quin ritme? de què depèn?

Pas 2:

:   $N_0=10$ cèules; $N$ és la variable depenent; $t$ serà la
    independent; n'hi ha alguna més? ($T$, aliment, posició a la placa,
    \...). Plantegem un model discret.

Pas 3:

:   Generem un diagrama de fluxe per entendre tot el que succeix en un
    intèrval de temps (`time step`)

    ::: center
    :::
:::

::: frame
(continued)

Pas 4:

:   Construïm les equacions:

    -   En $t=0$, $N=10$

    -   En $t=20$ minuts, $N=20=2\cdot N(0)$

    -   En $t=40$ minuts, $N=40=2 \cdot N(20)$, etc. $\Rightarrow b=2$

    $N(t=k\cdot t_1) = 2^k \cdot N_0$, amb $t_1=20$ minuts i
    $k\in\mathbf{N}$.

Pas 5:

:   Analitzem les equacions:

    :::: center
    ::: columns
       $t$   $N(t)$   $N(t)^{\mathrm{pred}}$
      ----- -------- ------------------------
       0.0    10.0             10.0
       1.0    18.8             20.0
       2.0    54.9             40.0
       3.0    94.1             80.0
       4.0   174.2            160.0

    ![image](ModelDiscretCultiuCellular.png){width="\\linewidth"}
    :::
    ::::
:::

::: frame
Exemple 1b: Model Continu Exemple d'un model continu per descriure
l'evolució del mateix cultiu ceular: Fórmula:
$N(t) = N_0 \cdot 2^{\frac{t}{20}}$ Comparació dels resultats dels
models continu i discret per diferents instants de temps.
:::

::: frame
Exemple 2: Model de Dolbear Exemple d'un model matemàtic que relaciona
la temperatura de l'aire ($T$) amb el nombre de carrisqueigs per minut
($N$) dels grills: Fórmula: $T = 50 + \frac{N-40}{4}$, on T és la
temperatura en F i $N$ el nombre de carrisqueigs per minut.
:::

# Limitacions

::: frame
Limitacions dels models

-   Els models matemàtics tenen limitacions, i els resultats són només
    tan interessants com les preguntes biològiques que els motiven.

-   Sovint, resoldre un model requereix fer suposicions dubtoses,
    fent-lo massa complicat.

-   Els models ens poden dir què és possible, però sense dades, no ens
    poden dir què ha passat o què està passant.

-   No té sentit promoure la biologia matemàtica per sobre d'altres
    àrees ni evitar completament les matemàtiques.

-   El progrés científic serà més ràpid gràcies a la unió entre la
    biologia matemàtica i l'empírica.
:::

# Conclusió

::: frame
Consell Final Practiqueu i feu exercicis! Igual que aprendre a jugar al
tennis, l'aprenentatge de les matemàtiques requereix pràctica. A classe
s'explicaran conceptes, però la pràctica és clau per assolir-los.
Afronteu l'assignatura amb predisposició i ganes per esdevenir
competents professionals en l'àmbit de la Biologia.
:::

::: frame
Bibliografia
:::
