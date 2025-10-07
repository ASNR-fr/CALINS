<PARTIE NECESSAIRE POUR LEXPORT DU README AU FORMAT PDF>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({ tex2jax: {inlineMath: [['$', '$']]}, messageStyle: "none" });
</script>

CALINS : CALculations and Investigations on Nuclear data uncertainties and Sensitivities

- [**INSTALLATION**](#installation)
- [**THÉORIE**](#théorie)
  - [**Introduction**](#introduction)
  - [**Vecteur de sensibilité**](#vecteur-de-sensibilité)
  - [**Matrice de covariances**](#matrice-de-covariances)
  - [**Règle d'alignement pour les calculs matriciels**](#règle-dalignement-pour-les-calculs-matriciels)
  - [**Calcul d'incertitude sur la réponse : σrespND**](#calcul-dincertitude-sur-la-réponse--σrespnd)
  - [**Assimilation de données expérimentales et calcul de biais par GLLSM**](#assimilation-de-données-expérimentales-et-calcul-de-biais-par-gllsm)
    - [**Pré-tri des cas benchmarks à assimiler**](#pré-tri-des-cas-benchmarks-à-assimiler)
    - [**Formules d'assimilation**](#formules-dassimilation)
  - [**Réflexions sur la théorie**](#réflexions-sur-la-théorie)
- [**MANUEL DEVELOPPEUR**](#manuel-developpeur)
  - [**methods.py**](#methodspy)
    - [**Acquisition des données utiles et formattage (Pandas Dataframe)**](#acquisition-des-données-utiles-et-formattage-pandas-dataframe)
    - [**Construction des matrices et vecteurs (Numpy Arrays)**](#construction-des-matrices-et-vecteurs-numpy-arrays)
    - [**Fonctions de calcul**](#fonctions-de-calcul)
    - [**Fonctions utilitaires**](#fonctions-utilitaires)
  - [**classes.py**](#classespy)
    - [**Case**](#case)
    - [**Uncertainty**](#uncertainty)
    - [**Bias**](#bias)
    - [**Assimilation**](#assimilation)
  - [**plots.py**](#plotspy)
  - [**errors.py**](#errorspy)


## **⚠️ Points de vigilance ⚠️**
(ce chapitre est à lire en connaissance du chapitre théorique et des commandes USER_EXAMPLES.ipynb)
- Certains isotopes particuliers présents dans les fichiers de sensibilités ne sont pas présents dans les matrices de variances-covariances, mais TSURFER (SCALE) utilise les données de variances-covariances d'autres isotopes pour les propager sur ces sensibilités, par exemple :
  - les isotopes liés (détectés par CALINS si la centaine de l'ID est >= 3), exemples : 
    - be-9 bounded (ID=4309 ou 4509) --> associé aux covariances de be-9 (ID=4009)
    - h-1 bounded (ID=1901 ou 1801) --> associé aux covariances de h-1 (ID=1001)
  - les isotopes avec une forme naturelle évaluée, exemple :
    - c-12 (bounded ?) (ID=6312) --> associé aux covariances de c-0 (ID=6000)
- Certaines matrices de covariances possèdent des données pour la réaction CAPTURE (ID=101), et d'autres pour la réaction N,GAMMA (ID=102) très proche de la catpure. Il faut donc être vigilant sur le choix des réactions à inclure dans les calculs (argument *reac_list*), soit 101, soit 102 en fonction de la matrice utilisée.
- Les fichiers *sdf* issus de calculs MCNP peuvent contenir des réactions avec un ID négatif :
  - CALINS associe la réaction -2 à la réaction 101, si le flag *mcnp* est activé dans la création de l'objet Case. Les autres réactions négative déclenchent simplement un warning.
- Si le fichier *sdf* possède plusieurs profils de sensibilité pour le même couple isotope-réaction (même IDs), TSURFER (SCALE) ne semble pas traiter *parser* les profils de la même manière entre SCALE-6.2.1 et SCALE-6.3.1. SCALE-6.2.1 semble faire une somme de tous les profils de sensibilité avec les même IDs isotope-réaction (ou bien prendre en compte seulement le dernier profil rencontré dans le fichier *sdf*) et SCALE-6.3.1 semble prendre en compte seulement le premier profil rencontré pour chaque couple isotope-réaction, et ignorer les autres si ils existent. Cela reste à vérifier. En revanche, CALINS possède 3 manières différentes de traiter ces occurrences de couple isotope-réaction, afin de pouvoir s'adapter au utilisateurs : 
  - *occurrences_rule = first* : la première occurrence du couple isotope-réaction (dans l'ordre de lecture du fichier) est utilisé comme profil de sensibilité pour cet isotope-réaction, et toutes les autres occurrences sont ignorées.
  - *occurrences_rule = sum* (option par défault): pour toutes les occurrences d'un couple isotope-réaction rencontrées, la somme des profils est utilisée comme profil de sensibilité pour cet isotope-réaction.
  - *occurrences_rule = first* : la dernière occurrence du couple isotope-réaction (dans l'ordre de lecture du fichier) est utilisé comme profil de sensibilité pour cet isotope-réaction, et toutes les autres occurrences sont ignorées.



# **INSTALLATION**

Après avoir télécharger le dossier github, placez-vous à sa racine (s'il n'est pas compressé), puis déclarez la variable d'environnement du proxy pour Pip : 

```
export PIP_PROXY=http://USER:PASSWORD@wpadirsn.proton.intra.irsn.fr:8088
```
L'utilisation du proxy est nécessaire seulement si vous êtes connectés au réseau IRSN. Le USER et PASSWORD doivent être en format URL (pas de charactères spéciaux de manière explicit, ou bien précédés de \ ), pour que le proxy fonctionne.
L'utilisation de l'option *pip install --proxy=...* ne fonctionne pas avec ce type de package, comprenant un fichier *pyproject.toml*, car l'option de proxy ne se transmet dans les subprocess d'installation de pip.
Si vous êtes sur Windows, utilisez la commande équivalente à *export ...* : *set ...*

Puis installez le package *calins* comme ceci en étant placé à la racine (s'il n'est pas compressé) :
```
pip install .
```
Sinon installez en vous plaçant dans le dossier contenant l'archive compressée du package :
```
pip install nom-de-larchive.zip
```

Vous pouvez vérifier que le package a bien été installé en vérifiant l'auto-completion dans un script python.


# **THÉORIE**

⚠️Le therme "réponse" (ou "resp" pour "response") sera employé ici pour désigner la réponse intégrale (expérimentale ou calculée) de manière générale⚠️


## **Introduction**
Quand une valeur de réponse est calculée, par exemple un k<sub>eff</sub>, il existe un biais de calcul Δresp (delta-réponse) entre la valeur de la réponse du modèle simulé et celle de l’objet physique. Ce biais est majoritairement constitué du biais dû aux données nucléaires Δ<sub>resp</sub><sup>ND</sup> et du biais dû aux approximations faites dans les schémas de calcul Δ<sub>resp</sub><sup>SC</sup>. Les méthodes décrites permettent d’estimer le biais de calcul de la réponse et son incertitude associée. 
La méthode d’estimation du biais de calcul Δ<sub>resp</sub> présentée ici est la GLLSM (méthode des moindres carrés linéaires généralisée). Elle repose sur la première hypothèse que le biais de calcul total est majoritairement constitué du biais Δ<sub>resp</sub><sup>ND</sup>. Cette hypothèse est vérifiée notamment dans le cas d’un calcul Monte Carlo ponctuel. Enfin, la seconde hypothèse est que la variation de la réponse engendrée par de faibles variations des données nucléaires est linéaire. 
La méthode de calcul de l’incertitude sur le calcul de la réponse σ<sub>resp</sub><sup>ND</sup> dû aux données nucléaires est un calcul matriciel nommé « Formule sandwich ». Il consiste en la propagation de sensibilités de la réponse aux données nucléaires pour un cas d’étude particulier, au travers d’un matrice d’incertitudes sur les données nucléaires (« matrice de covariances »).

A noter que l'on parle ici de "données nucléaires", mais que les méthodes présentées et les applications concernent seulement des données de sections efficaces microscopiques, à titre d'exemple de données nucléaires. On pourrait y inclure d'avantages de données comme la distribution angulaire de certaines réactions.

Enfin, il faut noter que les formules présentées ici sont valables pour des données exprimées en valeurs relatives.  
 

## **Vecteur de sensibilité**
Le vecteur de sensibilité, associé à une unique modélisation, est constitué d’une suite de sous-vecteurs de sensibilités de resp<sup>calc</sup> à une section efficace multi-groupes pour une réaction et un isotope particulier, comme décrit ci-dessous.
![](./images/Image_maker/Diapositive1.PNG)
Les paires isotope-réaction liées aux sous-vecteurs constituent une liste appelée « liste isotope-réaction ». Cette liste, et son ordre, sont importants pour le respect de la règle d’alignement des calculs matriciels décrits ici.

## **Matrice de covariances**
La matrice de covariances (Cov) est construite à partir des données d’incertitudes sur les valeurs des sections efficaces multi-groupes et des données de corrélation entre ces incertitudes contenues dans la matrice de corrélation (Corr). Par construction cette matrice est symétrique.
![](./images/Image_maker/Diapositive2.PNG)
![](./images/Image_maker/Diapositive3.PNG)
La matrice de covariances peut aussi être interprétée comme une composition de sous-matrices de covariances, pour chaque couple de paire (i<sub>h</sub>, r<sub>h</sub>)-(i<sub>v</sub>, r<sub>v</sub>), formés entre la liste isotope-réaction sur l’axe horizontal (i<sub>h</sub>, r<sub>h</sub>) et la liste isotope-réaction sur l’axe vertical (i<sub>v</sub>, r<sub>v</sub>). 
![](./images/Image_maker/Diapositive4.PNG)
## **Règle d'alignement pour les calculs matriciels**
Afin que chaque multiplication matricielle soit valide, il est nécessaire de respecter l'alignement des termes des deux matrices/vecteurs mis en jeux dans le calcul. Ainsi, l'application linéaire que constiue le produit matriciel doit multiplier des termes provenant du même couple isotope-réaction et du même groupe énergétique seulement. Un moyen simple de respecter ce critère, est de construire les vecteurs/matrices de sorte à ce qu'ils suivent la même liste isotope-réaction, dans le même ordre, avec des groupes d'énergies aussi dans le même ordre.

## **Calcul d'incertitude sur la réponse : σ<sub>resp</sub><sup>ND</sup>**
La formule sandwich permet la propagation de sensibilités de la réponse aux données nucléaires pour un cas d’étude particulier, au travers de la matrice de covariances, afin d’évaluer l’incertitude générale sur la réponse dû au données nucléaires σ<sub>resp</sub><sup>ND</sup> (en valeur relative).

**Formule Sandwich :**

$$ {σ_{resp}}^{ND} = \sqrt{S \space · \space Cov \space · \space S^t} $$

$S$ : vecteur de sensibilité &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$Cov$ : matrice de covariances

## **Assimilation de données expérimentales et calcul de biais par GLLSM**
La méthode GLLSM consiste en l'assimilation de données de simulation d'expérimentation (aussi appelés cas "benchmark") dont la réponse a été mesurée et calculée.

---
*Manuel SCALE-6.3.1* :
L'approche GLLS prend en compte les variations potentielles des données nucléaires et des réponses intégrales mesurées qui minimisent les différences entre les réponses intégrales mesurées et calculées (telles que keff) pour un ensemble d'expériences de référence, en tenant compte des incertitudes et des corrélations dans les données nucléaires évaluées et dans les mesures intégrales. Étant donné qu'il y a généralement beaucoup plus de valeurs de sections efficaces que de réponses intégrales mesurées, la détermination des modifications de données constitue un problème *under-determined*. Si les variations de données ne sont pas effectuées de manière raisonnable, des changements non physiques peuvent être obtenus. Les ajustements de données dans la méthodologie GLLS sont contraints par l'amplitude des incertitudes des données nucléaires et leurs corrélations. TSURFER/CALINS prend également en compte les incertitudes et les corrélations *(si fournies)* dans les mesures de la réponse intégrale, résultant des incertitudes expérimentales dans des paramètres tels que l'enrichissement, la densité, les impuretés, etc. Tant que des covariances de données réalistes sont utilisées, les ajustements de données résultants peuvent être considérés comme les « meilleures estimations » - dans les limites de l'approximation linéaire GLLS - pour des altérations réalistes des données qui améliorent les réponses intégrales calculées. Il peut être démontré que les équations GLLS fournissent l'estimateur du maximum de vraisemblance pour les données nucléaires correctes si les données nucléaires évaluées et les paramètres du système obéissent à une distribution normale multivariée [^1].
[^1]: R. N. Hwang. Topics in data adjustment theory and applications. In Proceedings of the Specialists' Meeting on the Application of Critical Experiments and Operating Data to Core Design via Formal Methods of Cross Section Data Adjustment. Jackson Hole, WY (USA), 9 1988.

(...)

La procédure GLLS peut modifier la valeur calculée pour l'application si celle-ci est « similaire » à certaines des réponses expérimentales. Dans ce cas, la réponse de l'application partage des caractéristiques de sensibilité aux données nucléaires similaires avec une ou plusieurs des réponses des expériences, et sera donc indirectement affectée par les mêmes ajustements de données qui influencent les réponses des expériences similaires. Cela fournit une méthode systématique et bien définie pour exploiter les mesures expérimentales de benchmark afin d'établir un biais et une incertitude dans le calcul de la réponse de l'application.

---

Cette méthode d'assimilation produit, d'une part, un vecteur de variations globales des données nucléaires (d'après les cas benchmarks assimilés) $\Delta\mu_{XS}$ permettant de prédire le biais de calcul de la réponse ("biais à posteriori") Δresp<sup>post</sup> pour un cas d'étude dont on ne connait pas la réponse de manière expérimentale (vecteur $S_{cas}$). D'autre part, l'assimilation produit une matrice de covariances ajustée ($Cov'$) permettant de calculer l'incertitude à posteriori σ<sub>resp</sub><sup>ND post</sup> d'un cas d'étude ($S_{cas}$) avec la formule sandwich. 

Le domaine de validité de cette méthode est délimité par plusieurs hypothèses : 
- le biais de calcul total est majoritairement constitué du biais Δ<sub>resp</sub><sup>ND</sup>. Cette hypothèse est vérifiée notamment dans le cas d’un calcul Monte Carlo ponctuel ;
- la variation de la réponse engendrée par de faibles variations des données nucléaires est linéaire.
- le cas d'étude est proche en terme de sensibilité des benchmarks assimilés (on parle de similarité)

### **Pré-tri des cas benchmarks à assimiler**
Les cas benchmarks assimilés ont d'autant plus d'impact sur la réduction de l'incertitude à posteriori σ<sub>resp</sub><sup>ND post</sup> du cas d'étude, qu'ils sont similaires à ce cas. Cela se traduit généralement par une proximité en terme de matière fissile, de type de modérateur, de spectre énergétique des neutrons, et du contexte physique en général.

Il est possible de faire un pré-tri des expériences à assimiler sur la base de ces critères physiques, i.e. les similarités entre les cas benchmarks et le cas d'étude. Par exemple, la classification des experiences de l'ICSBEP permet déjà de choisir parmis un lot d'expériences (cas benchmarks) similaire au cas d'étude en terme de matière fissile, ou de spectre.

On peut également définir des indicateurs de similarité entre deux vecteurs de sensibilités. Ils sont utiles pour faire un pé-tri des cas benchmarks qui, une fois assimilés, auront le meilleur impact sur l'incertitude à posteriori vis à vis du cas d'étude.
Ces indicateurs sont à calculer pour chaque benchmark que l'on souhaite ajouter à la liste de cas à assimiler. 

**Indice de similarité $E$ (valeur entre 0 et 1) :**
Cet indice est un produit scalaire normalisé entre un benchmark et le cas d'étude.

$$ E = {{S_1} \space · \space {S_2}^t \over ||S_1|| \space · \space ||S_2||} $$

**Indice de similarité $C_k$ (valeur entre 0 et 1) :**
Cet indice pondère le produit scalaire avec les valeurs de covariances.
$$ C_k = \sqrt{({{S_1} \space · \space Cov \space · \space {S_2}^t})^2 \over ({S_1} \space · \space Cov \space · \space {S_1}^t)({S_2} \space · \space Cov \space · \space {S_2}^t)} $$

**Indice de recouvrement G (valeur entre 0 et 1) :**
Cet indice calcul le taux de recouvrement de chaque sensibilité correspoondant à un groupe d'énergie, un isotope et une réaction. Ce taux est calculé en prenant un cas d'étude en tant que cas de référence, et l'autre cas d'étude comme cas comparé.

$$  SS = 1 - {\sum_{g}^{}\sum_{i, r}^{} \left\{\begin{matrix}
if \space \space {S_{ref}^{g, i, r}}\times {S_{comp}^{g, i, r}} > 0 \space and \space \left|{S_{ref}^{g, i, r}}\right|\geq \left|{S_{comp}^{g, i, r}}\right| : {S_{ref}^{g, i, r}} - {S_{comp}^{g, i, r}} \\
else \space : {S_{ref}^{g, i, r}}
\end{matrix}\right. \over 
\sum_{g}^{}\sum_{i, r}^{} {S_{ref}^{g, i, r}}  }$$

${S_{ref}}$ : vecteur de sensibilité de référence &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ${S_{comp}}$ : vecteur de sensibilité comparé


**Indice de recouvrement SS - Shared Sensitivity (formule de Mariya BROVCHENKO) :**
Cet indice calcul le taux de recouvrement de chaque sensibilité correspoondant à un groupe d'énergie, un isotope et une réaction. Ce taux est calculé en prenant un cas d'étude en tant que cas de référence, et l'autre cas d'étude comme cas comparé.

$$ SS =  {\sum_{g}^{}\sum_{i, r}^{} \left\{\begin{matrix}
if \space \space {S_{ref}^{g, i, r}}\times {S_{comp}^{g, i, r}} > 0 \space : min(\left|{S_{ref}^{g, i, r}}\right|, \left|{S_{comp}^{g, i, r}}\right|) \\
else \space : 0
\end{matrix}\right. \over 
\sum_{g}^{}\sum_{i, r}^{} \left| {S_{ref}^{g, i, r}} \right| }$$

${S_{ref}}$ : vecteur de sensibilité de référence &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ${S_{comp}}$ : vecteur de sensibilité comparé

**Tri sur la valeur du $\chi ^2$**

*A rédiger : cf. Chaptire 6.8.4.2.1 du Manuel Utilisateur SCALE 6.3.1 [^2]*
[^2]: https://scale-manual.ornl.gov/tsurfer.html#consistency-relations-and-chi-square-filtering

### **Formules d'assimilation**
**Description des données utilisées pour la GLLSM** : 
- un "vecteur de C-E benchmarks" $\Delta k_{C-E \space bench}$ qui contient les écarts relatifs sur la réponse ((C-E)/C) des benchmarks retenus ;
![](./images/Image_maker/Diapositive5.PNG)
- les incertitudes σ<sub>resp</sub><sup>exp</sup>, pour construire la "matrice d'incertitudes expérimentales" $C_{bench}$. Cette matrice pourrait intégrer les corrélations experimentales entre les différents benchmarks. Ce n'est pas le cas ici car les données manquent pour réaliser cela (ce point est aborder plus bas dans le chapitre "Réflexions sur les corrélations expérimentales"). ;
![](./images/Image_maker/Diapositive6.PNG)
- la "matrice des sensibilités benchmarks" $S_{bench} = ..|S_{bench}^j | S_{bench}^{j+1} | ..$ construite à partir des vecteur de sensibilité $S_{bench}^j$ (j : indice du benchmark concerné)
![](./images/Image_maker/Diapositive7.PNG)

⚠️Il est à noter que la matrice $Cov$ dépend des benchmarks utilisés. En effet elle est définies par l'union des couples isotope-réaction présent dans les benchmarks et le cas d'étude (point précisé dans le chapitre "Manuel développeur").⚠️
***
**Vecteur de variations globales des données nucléaires (d'après les cas benchmarks assimilés) $\Delta\mu_{XS}$ :** 

$$ \Delta\mu_{XS} = - Cov \space · \space {S_{bench}}^t \space · \space  (C_{bench}+S_{bench} \space · \space Cov \space · \space {S_{bench}}^t)^{-1} \space · \space \Delta k_{C/E ~ bench} $$


**Matrice de covariances ajusté $Cov'$ :**

$$ Cov' = Cov - \Delta Cov_{assim} $$

$$ Cov' = Cov - Cov \space · \space {S_{bench}}^t \space · \space  (C_{bench}+S_{bench} \space · \space Cov \space · \space {S_{bench}}^t)^{-1} \space · \space  S_{bench} \space · \space Cov $$ 


L'inversion necessaire n'est pas toujours mathématiquement possible. L'impact des méthodes d'inversion conduisant à des termes "pseudo-inverse" n'a pas encore été étudié. Pour vérifier l'inversibilité, la méthode algorithmique "Cholesky" d'inversion de matrice conduit soit à un résultat valide, soit à un arrêt de l'algorithme. 
***
**Biais à posteriori Δresp<sup>post</sup>**, écart $(k_{eff}^{expe}-k_{eff}^{calc}) \over k_{eff}^{calc}$ (unité : %, valeur relative à la réponse calculée) :

$$ \Delta {resp}^{post} = {S_{cas}}  \space · \space  \Delta\mu_{XS} $$


**Incertitude à priori** (unité : %, valeur relative à la réponse calculée) :

$$ {σ_{resp}}^{ND~ prior} = \sqrt{{S_{cas}} \space · \space Cov \space · \space {S_{cas}}^t} $$

**Incertitude à posteriori** (unité : %, valeur relative à la réponse calculée) :

$$ {σ_{resp}}^{ND~ post} = \sqrt{{S_{cas}} \space · \space Cov' \space · \space {S_{cas}}^t} $$

***

## **Réflexions sur les corrélations expérimentales**
Compte tenue de la méthodologie présenté ci-dessus, le fait de fixer les corrélations expérimentales à 0 a plusieurs conséquences. Il a été observé qu'une assimilation comportant plusieurs fois un cas benchmarks exactement identique, obtenait une incertitude a posteriori plus faible, et un  biais différent qu'une assimilation composée d'une seule fois ce cas benchmark. Cet effet est dû au postulat implémenté que ces multiples expériences ont une corrélation nulle. En effet, le fait qu'une même expérience a réellement eu lieu plusieurs fois, et que toute ses réalisations soient assimilées, ne devraient pas avoir d'impact sur le calcul d'incertitude du cas d'étude, car aucune information n'a réellement été ajouté mis à part le fait que cette expériences à été réalisée plusieurs fois - donnée qui ne devrait pas avoir d'impact sur notre cas d'étude, sauf dans le cas où l'on souhaite attribuer plus de valeurs à une expérience qui a été plus reproduite qu'une autre. 

 Cela questionne sur la ressemblance entre les cas benchmarks assimilés quand la méthode GLLS est employée. Si plusieurs expériences ont intrinsèquement des sensibilités, et des réponses (expérimentales et calculées) extrêmement similaires (expériences d'une même série par exemple) alors le non-ajustement des corrélations biaise le processus d'assimilation. Cela semble avoir pour effet de sur-estimer l'incertitude a posteriori et de mal-évaluer le biais correspondant.

 Des études sur l'impact des corrélations expérimentales ont déjà été menée notamment [^3].
 [^3]: T. Nicol, C. Carmouze. Impact of experimental correlation on transposition method carry out with critical integral experiments. ICNC 2019 - 11th International conference on Nuclear Criticality Safety,Sep 2019, Paris, France. ffcea-02614125


# **MANUEL DEVELOPPEUR**
Les multiplications matricielles sont effectuées avec l'opérateur de python "A @ B" étant totalement équivalent à l'opérateur numpy.matmul.

⚠️Attention : la priorité opératoire est trés importante quand il s'agit des méthodes de mutliplications matricielles algorithmiques, afin d'optimiser le temps de calcul.⚠️

---
## **methods.py**
### **Acquisition des données utiles et formattage (Pandas Dataframe)**
Fonctions concernées : *format_ ... _to_dataframe (...)*

Ces fonctions ont pour but de lire les différents fichiers contenant les matrices de covariances et les sensibilités, afin de formatter ces données dans des DataFrame Pandas, faciles à manipuler pour construire les matrices et vecteurs.

Les DataFrame de **sensibilités** sont construit seulement à partir de fichiers au format *.sdf* comme ceci : 

| ISO | REAC | SENSI |
| :-: | :-: | :-: |
|iso1|reac1[iso1] | $S^{1, 1}$
|iso1|reac2[iso1] | $S^{1, 2}$
|iso2|reac1[iso2] | $S^{2, 1}$
|iso2|reac2[iso2] | $S^{2, 2}$
|.|. | .
|.|. | .

Avec :
- iso i : numéro d'identification de l'isotope i (Z x 1000 + A), au format *string* ou *int* ;
- reac r[iso i] : numéro d'identification de la réaction r (en fonction du dictionnaire *reac_trad*), indépendant de l'itosope i ;
- $S^{i, r}$ : liste python [. , . , ..] comprenant les sensibilités relatives multi-groupes de la réponse à la réaction r de l'itosope i, en ordre décroissant d'énergies ;

***Remarque*** : tous les isotopes et réactions rencontrés sont stocké dans ce dataframe ; dont la réaction TOTAL - ID : 1. Cependant un filtre sur la réaction TOTAL est appliqué lors de la construction des vecteurs *numpy* de sensibilité, avant leur opérations.

Les DataFrame des **matrices de covariances** sont construit comme ceci : 

| ISO_H | REAC_H | ISO_V | REAC_V | STD |
| :-: | :-: | :-: | :-: | :-: |
|iso1_h|reac1_h[iso1_h] | iso1_v|reac1_v[iso1_v] |$Cov^{(1, 1), (1, 1)}$
|iso1_h|reac2_h[iso1_h] | iso1_v|reac2_v[iso1_v] |$Cov^{(1, 2), (1, 2)}$
|iso1_h|reac1_h[iso1_h] | iso1_v|reac2_v[iso1_v] |$Cov^{(1, 1), (1, 2)}$
|iso2_h|reac2_h[iso2_h] | iso2_v|reac2_v[iso2_v] |$Cov^{(2, 2), (2, 2)}$
|.|.|.|.|.
|.|.|.|.|.

Avec :
- iso i_h : numéro d'identification de l'isotope i_h (Z x 1000 + A), au format *string* ou *int*, dont les groupes d'énergies sont décroissants sur l'axe **horizontal** de la sous-matrice $Cov^{(i_h, r_h), (i_v, r_v)}$ ;
- reac r_h[iso i_h] : numéro d'identification de la réaction r_h (en fonction du dictionnaire *reac_trad*), dont les groupes d'énergies sont décroissants sur l'axe **horizontal** de la sous-matrice $Cov^{(i_h, r_h), (i_v, r_v)}$ ;
- iso i_v : numéro d'identification de l'isotope i_v (Z*1000 + A), au format *string* ou *int*, dont les groupes d'énergies sont décroissants sur l'axe **vertical** de la sous-matrice $Cov^{(i_h, r_h), (i_v, r_v)}$ ;
- reac r_v[iso i_v] : numéro d'identification de la réaction r_v (en fonction du dictionnaire *reac_trad*), dont les groupes d'énergies sont décroissants sur l'axe **vertical** de la sous-matrice $Cov^{(i_h, r_h), (i_v, r_v)}$ ;
- $Cov^{(i_h, r_h), (i_v, r_v)}$ : liste 2D (liste de liste) python [[. , . , ..], [. , . , ..], ..] contenant la sous-matrice de covariances multi-groupes entre la réaction r_h de l'isotope i_h et la réaction r_v de l'isotope i_v ;

Il existe, entre autre, une fonction format_..._to_dataframe(...) pour chaque matrice de covariances d'origine différente. Les matrices lisisbles pour l'instant sont :
- SCALE au format AMPX (sous la forme d'un fichier)
- SCALE au format binaire (sous la forme d'un fichier)
- COMAC (sous la frome d'un dossier de fichiers)
- GENDF (*work in progress*)

**ATTENTION** : La **norme** de *CALINS* est de formatter les données de sensibilités et de covariances dans l'ordre **DECROISSANT** de groupes d'énergies sur la gamme 20 MeV -> 0 MeV.

### **Construction des matrices et vecteurs (Numpy Arrays)**
Fonctions concernées : *make_ ... (...)*

Ces fonctions permettent de construire les vecteurs et matrices sous forme de *Numpy Array* en tenant compte de la **règle d'alignement** (cf. chapitre *Règle d'alignement pour les calculs matriciels*), afin que les vecteurs/matrices soient prêts pour les operations matricielles.

Pour que la règle d'alignement soit respectée entre plusieurs vecteurs de sensibilités et plusieurs matrices de covariances, ils sont construit simultanément sur la base d'une liste isotope-réaction commune. Cette liste dictera l'ordre des termes des vecteurs de sensibilités, ainsi que le placement des sous-matrices dans la matrice complète de covariances.
Cette liste se construit d'abord comme l'intersection ou l'union des listes isotope-réaction propres à chaque vecteur de sensibilité mis en jeu (en fonction des operations futurs que vont subir ces vecteurs - par défault l'option 'union' est retenue) - fonction *get_common_iso_reac_list* ; puis si une matrice de covariance est aussi mise en jeux, la liste isotope-reaction est réduite à son intersection avec la liste des isotopes-réactions présents dans les données de la matrice de variances-covariances (cette opération est effectué dans la fonction *make_cov_matrix*). 

**ATTENTION** : à savoir que si des isotopes non-conventionnels comme le "h-poly" ID=1901 est présent dans la liste isotope-réaction issue des vecteurs de sensibilités, et non dans le DataFrame de la matrice de variances-covariances, alors le DataFrame est enrichi, de sorte à avoir des données de variances-covariances pour cette isotope, si son isotope de base (ici h-1 ID=1001) est déjà présent dans le DataFrame. Les données utilisées pour cette isotope non-conventionnel sont une copie des données déjà présente pour l'isotope de base. La méthode actuelle pour détecter des isotopes non-conventionnels est de vérifié si la centaine de son ID est supérieur ou égale à 3 (exemple : h-poly / 1901, 9>3), d'après nos observations.
Le même principe est appliqué pour les isotopes possédant une forme "naturel" (ID=xx000) dans les données de variances-covariances. Si un isotope issue des vecteurs de sensibilités n'est pas présent dans les données du DataFrame de variances-covariances, alors on vérifie si des données pour sa forme naturel sont déjà présentes. Si oui, le DataFrame est enrichi avec une copie de ces données, associées à l'isotope (et ses réaction) non-présents initialement.

***Remarque :*** Un filtre excluant la réaction TOTAL est appliqué à cet endroit sur la liste isotope-réaction.

Pour les vecteurs de sensibilités, cette liste est parcourue dans l'ordre, et pour chaque paire isotope-réaction l'on recherche si des données concernant cette paire sont présentes dans son Dataframe. Si oui, alors les sensibilités multi-groupes sont ajoutée à la suite de son Numpy Array en construction. Si non, alors des sensibilité **nulles** multi-groupes sont ajoutée à la place.

Pour la matrice, elle est initialisée par une matrice de zéros, de dimmensions 2D *[(N_group x Nb_iso-reac),  (N_group x Nb_iso-reac)]*. Puis, la liste isotope-réaction est parcourue deux fois dans deux boucles *for* imbriquées (pour les indices horizontal et vertical). Et la sous-matrice du DataFrame correspondant aux deux paires isotope-réaction (i_h, r_h) et (i_v, r_v), ($Cov^{(i_h, r_h), (i_v, r_v)}$), est insérée dans la matrice $Cov$ aux bonnes coordonnées, si elle est bien présente dans le DataFrame. On vérifie aussi si la sous-matrice inverse correspondant aux deux paires (i_h, r_h) et (i_v, r_v) **interchangées** ($Cov^{(i_v, r_v), (i_h, r_h)}$) existe dans son DataFrame, et dés lors la **transposée** de la sous-matrice est insérée dans la matrice aux bonnes coordonnées. Enfin, si aucune des deux correspondances existe dans le DataFrame, alors la sous-matrice reste nulle à ces coordonnées. 

### **Fonctions de calcul**
Fonctions de type : *calcul_... (...)* 

Ces fonctions utilisents les fonctions *make_...* afin de construire les vecteurs de sensibilités et la matrice de covariance sur la base d'une même liste isotope-réaction, pour faire les calculs matriciels.

Il existe une fonction par indice de similarité. Les indices implémentés sont les suivant :
- E
- C<sub>k</sub>
- indice de recouvrement SS (formule du CEA)
- indice de recouvrement SS (formule de Mariya BROVCHENKO)

Ces fonctions ont pour input au moins deux vecteurs de sensibilité et une matrice de covariances dans le cas du calcul du C<sub>k</sub>. Les vecteurs de sensibilité peuvent être soit des chemins vers des fichiers .sdf, soit des objet *case* (présenté plus loin), mais pas de vecteurs Numpy déjà construits.

Enfin une troisième fonction permet le calcul d'une incertitude, avec en input un vecteur de sensibilité et une matrice de covariances. Le vecteur de sensibilité peut être un chemin vers un fichier *.sdf*, ou un objet *case*. Cette fonction ne renvoie pas explicitement un résultat mais renvoie l'objet *uncertainty*. Cette objet est décrit dans le chapitre ***classes***.

### **Fonctions utilitaires**
Le module *methods* contient aussi des fonctions utilitaires qui permettent de convertir un numéro d'identification isotopique en chaîne de charactère, et vis et versa.

## **classes.py**

### **Case**
L'objet *case* permet de construire un cas d'étude (ou cas benchmark), à partir d'un fichier de sensibilité *.sdf*, et de stocker les données utiles autour du cas. L'utilisation d'un objet rend donc ces informations accessibles facilement. Les propriétés disponibles sont les suivantes :
- *case.sdf_path* : le chemin vers le cas ;
- *case.casename* : le nom du cas d'étude (initialisé comme étant la base du chemin vers le fichier *.sdf*) ;
- *case.group_nb* : nombre de groupes d'énergies du maillage utilisé ;
- *case.e_bins* : bornes des groupes d'énergies du maillage utilisé ;
- *case.iso_reac_list* : liste isotope-réaction dont les données sont disponibles dans le fichier de sensibilité, sous la forme des numéros d'identifiant des isotopes et réactions ;
- *case.sensitivities* : DataFrame contenant les sensibilités du cas, obtenu grace à la fonction *format_sensi_to_dataframe(...)* (*methods.py*)
- *case.resp_calc* : valeur de la réponse calculé/modélisé ;
- *case.sigma_resp_calc* : incertitude sur la réponse calculée/modélisée ;
- *case.resp_expe* : valeur de la réponse expérimentale (pour les cas benchmarks) ;
- *case.sigma_resp_expe* :incertitude sur la réponse expérimentale (pour les cas benchmarks);

Enfin, l'objet *case* possède une fonction permettant d'afficher et d'enregistrer au format *.html* un histogramme des sensbilités intégrales, associé aux profils de sensibilité, pour chaque isotope et réaction isotopique.


### **Uncertainty**
L'objet *Uncertainty* permet de calculer une incertitude (absolue) à partir de la formule Sandwich, et de stocker plusieurs propriétés autour de ce calcul. Cet objet a pour inputs un vecteur de sensibilité et une matrice de covariances déjà construits sous format Numpy Array, ainsi que le nombre de groupe d'énergies et la liste isotope-réaction qui décrivent la construction du vecteur et de la matrice.  Les propriétés disponibles sont les suivantes :
- *Uncertainty.resp_calc* : valeur de la réponse calculée du cas d'étude ;
- *Uncertainty.value* : valeur de l'incertitude absolue calculée, en pcm ;
- *Uncertainty.group_nb* : nombre de groupes d'énergies du maillage utilisé ;
- *Uncertainty.e_bins* : bornes des groupes d'énergies du maillage utilisé ;
- *Uncertainty.iso_reac_list* : liste isotope-réaction complète, dans l'ordre de construction du vecteur et de la matrice, sous la forme des numéros d'identifiant des isotopes et réactions ;
- *Uncertainty.decomposition* : contribution à l'incertitude relative (%) au carré, pour chaque couple isotope-réaction, sous la forme d'un Dataframe. Cette décomposition correspond à la décomposition de la combinaison linéaire de la dernière étape du calcul matriciel de la formule Sandwich au carrée (X · S<sub>cas</sub>) comme décrit ci-dessous (elle prend donc en compte les co-variances aussi) :
$$({σ_{resp}}^{ND})^2 = {S_{cas}} \space · \space Cov \space · \space {S_{cas}}^t = X \space · \space {S_{cas}}^t$$
avec : $\space X =  {S_{cas}} \space · \space Cov \space \space$ (*DIM* : 1D -> I<sup>r</sup> x G)

Decompositon par iso-réac :  
$d^1 = (X^1 \space · \space S_{cas}^1) \space$ |
$\space d^2=(X^2 \space · \space S_{cas}^2) \space$ |
$\space d^3 = (X^3 \space · \space S_{cas}^3) \space $ |
$\space ...$

avec : 

$d^i = X^i \space · \space S_{cas}^i$

$S_{cas}^i$ : sous-vecteur de sensibilité pour le couple iso-réac *i*.

L'objet *Uncertainty* possède une fonctions permettant d'afficher et/ou enregistrer au format *.html* plusieurs éléments d'analyse : 
 - *Uncertainty.export_results()* :
   - l'incertitude ;
   - la valeur réponse calculée et son incertitude de schéma de calcul ;
   - le nombre de groupes d'énergies ;
   - le vecteur de sensibilité du cas d'étude ;
   - les décompositions par couple iso-réac de l'incertitude ;
   - les histogrammes des intégrales de la matrice de covariances par couple isotope-réacton ;
   - les sous-matrices de covariances pour des couples isotope-réaction choisis pas l'utilisateur ;

Enfin, l'objet *Uncertainty* est invoqué dans la fonction *calcul_uncertainty(...)*, et deux fois en tant que propriété de l'objet *Assimilation* (pour l'incertitude à priori et l'incertitude à posteriori).

### **Bias**
L'objet *Bias* est similaire à l'objet *uncertainty*. Il permet de calculer le biais à posteriori Δresp<sup>post</sup> à partir de la formule du biais (cf. chapitre *THEORIE*), et de stocker les mêmes propriétés que l'objet *uncertainty* à propos du calcul. Cet objet a pour input un vecteur de sensibilité et le vecteur de variations globales des données nucléaires $\Delta\mu_{XS}$ (après assimilation) déjà construits sous format Numpy Array, ainsi que le nombre de groupe d'énergies et la liste isotope-réaction qui décrivent la construction des vecteurs. Les propriétés disponibles de l'objet  sont les suivantes :
- *Bias.resp_calc* : valeur de la réponse calculée du cas d'étude ;
- *Bias.value* : valeur du biais absolu calculée, en pcm ;
- *Bias.group_nb* : nombre de groupes d'énergies du maillage utilisé ;
- *Bias.e_bins* : bornes des groupes d'énergies du maillage utilisé ;
- *Bias.iso_reac_list* : liste isotope-récation complète, dans l'ordre de constructon du vecteur et de la matrice, sous la forme des numéros d'identifiant des isotopes et réactions ;
- *Bias.decomposition* : contribution au biais relatif (%), pour chaque couple isotope-réaction (sous la forme d'un Dataframe) comme décrit ci-dessous :

$$ \Delta {resp}^{post} = {S_{cas}}  \space · \space  \Delta\mu_{XS} $$

Décomposition par iso-réac : $d^1 = S_{cas}^1 \times \Delta\mu_{XS}^1 \space ; \space d^2 = S_{cas}^2 \times \Delta\mu_{XS}^2 \space ; \space d^3 = S_{cas}^3 \times \Delta\mu_{XS}^3 \space ; \space ...$

avec : 

$d^i = S_{cas}^i  \space · \space \Delta\mu_{XS}^i$ 

$S_{cas}^i$ : sous-vecteur de sensibilité pour le couple iso-réac *i*.

Enfin, l'objet *Bias* est seulement invoqué en tant que propriété de l'objet *Assimilation*.

### **Assimilation**
L'objet *Assimilation* permet de réaliser la GLLSM d'après une liste de cas benchmarks, un cas d'étude et une matrice de covariances. Les cas benchmarks et le cas d'étude peuvent être soit sous forme de chemins vers les fichiers *.sdf*, soit sous la forme d'un objet *Case*. A contrario, la matrice de covariances doit être formaté en DataFrame (pour l'instant). Cette classe utilise les fonctions de *methods*, *plots*, et *errors* afin de formatter les données, construire les objets matrices et vecteurs, effectuer les calculs matriciels, et afficher les données utiles sous forme de plots.

L'initialisation d'un objet *Assimilation* est constitué de plusieurs étapes principales :
1. le formattage des données de sensibilités (benchmarks, et cas d'étude) en Dataframe ;
2. la construction de la matrice des sensibilités benchmarks, du vecteur de sensibilité du cas d'étude, de la matrice de covariances, suivant le même ordre isotope-réaction ;
3. la construction de la matrice d'incertitudes des réponses des benchmarks, et du vecteur de C-E benchmarks, suivant le même ordre des cas benchmarks ;
4. *self.check_dimensions()* & *check_correspondences()*: vérification de la conformité des dimensions de tous les vecteurs/matrices, et création de warning pour les sensibilités élevées qui n'ont pas de données de variances-covariances ;
5. *self.calcul_prior_uncertainty()* : calcul et création d'un objet *uncertainty* ;
6. *self.calcul_matrix_assimilation()* : calcul de $\Delta \mu_{XS}$ , et $\Delta Cov_{assim}$ et application du filtre $\chi ^2$ si l'utilisateur à donner une valeur cible, et application du filtre sur le Ck si l'utilisater à donner une valeur cible ;
7. *self.calcul_bias()* : calcul et création d'un objet *bias* ;
8. *self.calcul_post_uncertainty()* : calcul et création d'un objet *uncertainty* ;

L'objet *Assimilation* possède deux fonctions permettant d'afficher et/ou enregistrer au format *.html* plusieurs éléments d'analyse : 
 - le vecteur de sensibilité du cas d'étude
 - un fichier de sortie comprenant tous les paramètres de l'assimilation : 
   - les incertitudes à priori et à posteriori, ainsi que le biais ;
   - la valeur du $\chi ^2$ intial et final ;
   - le nombre de cas benchmarks retirés ;
   - le nombre de groupes d'énergies ;
   - la liste des cas benchmarks avec leurs chemins, les valeurs de la réponse calculées et expérimentales, les $\chi ^2_{individuel}$, les coefficients de similarité Ck avec le cas d'étude,et un fond rouge pour les cas exclus de l'assimilation ;
   - le vecteur de sensibilité du cas d'étude ;
   - les décompositions par couple iso-réac des deux incertitudes et du biais ;
   - les histogrammes des intégrales de la matrice de covariances par couple isotope-réacton, et de la matrice $\Delta Cov_{assim}$ en guise de comparaison ;
   - les sous-matrices de covariances et de $\Delta Cov_{assim}$ pour des couples isotope-réaction choisis pas l'utilisateur ;

## **plots.py**
Les fonctions présentes utilisent le package de visualisation graphique *Plotly* (v5.14), et notament son module plotly.express. Ce package est ergonomique et explicite, et permet de manipuler chaque graphique comme un lot de données html, afin de pouvoir combiner plusieurs  graphique dans un même fichier *.html* si besoin.

Les fonctions de *plots* sont destinées à une utilisation au sein des objets, mais peuvent-aussi être utilisées indépendamment par les utilisateurs. Plusieurs fonctions sont utilisées pour des types de graphiques différents : 
- un histogramme des sensbilités d'un cas d'étude, pour chaque isotope et réaction isotopique ; 
- un histogramme des covariances intégrales d'une matrice pour chaque isotope et réaction isotopique ;
- une *colormap* 2D d'une sous-matrice de covariances, pour une paire de couples (iso-réac)<sup>h</sup> - (iso-réac)<sup>v</sup> donnée ;
- un tableau au format *.html* d'après des données sous forme de listes ;

## **errors.py**
Les erreur customisées implémentées sont les suivantes :
- *EmptyParsingError* : appelée dans *methods*, dans le cas d'une acquisition et d'un fromattage de données qui aboutit à un/une vecteur/matrice vide (erreur de formattage/construction);
- *DimError* : appelée quand les dimensions de plusieurs vecteurs/matrices ne sont pas cohérentes entre elles pour des mutliplications matricielles ;
- *SensInputError* : quand des données de sensibilités n'ont pas le bon format pour l'acquisition et le formattage ; 

## **logs.py**
Les fonctions permettent de tracer les actions réalisée par l'utilisateur dans une optique de vérification et non dans une optique développeur de débogage. Les logs sont à la fois écrits sur le terminal (sauf si la variable interne logs.VERBOSE=False, alors le terminal ne recoit que les warnings) et dans un fichier de log situé dans *~USER/.CALINS/*. Ce dossier est checké à chaque import de CALINS, et vidé de son dernier fichier si le nombre de fichier dépasse 10 (car un nouveau fichier se créer chaque jour différent ou CALINS est importé).

Le décorateur *log_exec* permet de suivre l'exécution des fonctions principales de CALINS.   
