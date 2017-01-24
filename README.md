# Installation

Démarrez Rstudio, puis installez les paquets *bnlearn* et *Rgraphviz*.

```R
install.packages("bnlearn")
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "Rgraphviz"))
```

# Prise en main de *R*

Créez un nouveau projet: File -> New Project -> New Directory

Créez un nouveau script R dans lequel vous allez travailler: File -> New File -> R Script

Commencer par charger la base **alarm**:

```R
library("bnlearn")

# charge la base de données alarm
data("alarm")

# infos sur les colonnes (nos variables)
ncol(alarm)
colnames(alarm)

# infos sur les lignes (nos observations)
nrow(alarm)
rownames(alarm)

# dix premières lignes, 5 premières colonnes
alarm[1:10, 1:5]

# lignes 3, 5, 1, colonnes "ANES", "HIST" et "MINV"
alarm[c(3, 5, 1), c("ANES", "HIST", "MINV")]
```

Quelques commandes utiles:

+ Auto-complétion: *ctrl+espace*;
+ Page d'aide d'une fonction: `?nom_fonction` (à exécuter dans l'interpréteur);
+ Exécuter la sélection ou la ligne courante: *ctrl+entrée*;
+ Exécuter le fichier courant: *ctrl+shift+S*.

# Prise en main de *bnlearn*

Le package *bnlearn* permet, entre autres, de faire des tests d'indépendance conditionnelle.

```R
ci.test(x = "PAP", y = "SHNT", z = as.character(NULL), data = alarm, test = "mi")

res = ci.test(x = "PAP", y = "SHNT", z = "PMB", data = alarm, test = "mi")
res$statistic
res$p.value
```

Inspectez la relation entre les variables *PAP* et *SHNT*:

```R
table(alarm[, "PAP"])
plot(alarm[, "PAP"])
prop.table(table(alarm[, "PAP"]))

table(alarm[, "SHNT"])
plot(alarm[, "SHNT"])
prop.table(table(alarm[, "SHNT"]))

ct = table(alarm[, c("PAP", "SHNT")])
prop.table(ct)
prop.table(ct, margin = 1)
prop.table(ct, margin = 2)
```

Parmis les relations d'indépendance ci-dessous, lesquelles sont vraies?

+ `STKV ⟂ HR | ∅`;
+ `STKV ⟂ HR | CO`;
+ `HR ⟂ CO | ∅`;
+ `HR ⟂ CO | STKV`;
+ `CO ⟂ STKV | ∅`;
+ `CO ⟂ STKV | HR`.

Quelle structure de réseau Bayésien permet d'encoder le modèle d'indépendance entre les trois variables *STKV*, *HR* et *CO* ?

Inspectez la relation entre *STKV* et *HR*:

```R
mask = rep(TRUE, nrow(alarm))
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x)")
```

Inspectez la relation entre *STKV* et *HR* sachant *CO* (vous pouvez remplacer `"HIGH"` par `"LOW"` ou `"NORMAL"`):

```R
mask = alarm[, "CO"] == "HIGH"
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x,z=HIGH)")
```

# Inférence dans un réseau Bayésien

Tout d'abord, construisez un réseau Bayésien complet (structure et paramètres) avec les instructions suivantes:

```R
# structure
bn = hc(alarm)
graphviz.plot(bn)

# parametres
bn = bn.fit(bn, data = alarm, method = "bayes")
bn[["CO"]]
```

## Inférence approchée

On peut calculer (inférer) n'importe quelle probabilité à partir du réseau bayésien avec la commande `cpquery()`. Exécutez plusieurs fois les instructions suivantes. Qu'observez-vous?

```R
cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW"))
cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW" & CO == "LOW"))
```

## Inférence exacte

Récupérez le fichier [**includes.R**](https://github.com/gasse/pgm_tutorial/blob/master/includes.R) qui contient la fonction `exact.dist()`, et ajoutez-le à votre projet. Vous pouvez désormais faire de l'inférence exacte comme suit:

```R
source("includes.R")

p = exact.dist(bn, event = c("STKV", "HR", "CO"), evidence = TRUE)
sum(p["HIGH", "LOW", ]) / sum(p[, "LOW", ])
sum(p["HIGH", "LOW", "LOW"]) / sum(p[, "LOW", "LOW"])
```

Attention, l'inférence exacte dans un réseau bayésien peut être très coûteuse:

```R
p = exact.dist(bn, event = c("INT", "APL"), evidence = TRUE)
```

## Do-calculus

A l'aide de la fonction `exact.dist()`, calculez la distribution conditionnelle de *HYP* sachant *STKV*, autrement dit p(y|x). Puis, en supposant que le réseau bayésien est causal, calculez la distribution de probabilité de *HYP* sachant que la valeur de *STKV* a été forcée, autrement dit p(y|do(x)). On rappelle la formule d'ajustement pour "supprimer" une variable confondante: p(y|do(x)) = ∑<sub>z</sub> p(y|x,z)p(z)

# L'algorithme PC

Nous allons maintenant coder une algorithme d'apprentissage de structure simple: l'algorithme `PC`.

## Construire le squelette

Commencez par initialiser un squelette (graphe non-dirigé) complet. A partir d'un graphe vide, ajoutez un arc non-dirigé entre chaque paire de noeuds:

```R
vars = colnames(alarm)
g = empty.graph(vars)
for (x in vars) {
  for (y in setdiff(vars, x)) {
    g = set.edge(g, from = x, to = y)
  }
}
graphviz.plot(g)
```

Ensuite, épurez le squelette en suivant l'algorithme suivant:

1. pour chaque paire de variables X et Y, tester X ⟂ Y | ∅. Si la relation est vraie, alors retirer l'arc correspondant;
2. pour chaque paire (ordonnée) de variables X et Y, et pour chaque variable Z adjacente à X, tester X ⟂ Y | Z. Si la relation est vraie, alors retirer l'arc correspondant;
3. pour chaque paire (ordonnée) de variables X et Y, et pour chaque ensemble **Z** de 2 variables adjacentes à X, tester X ⟂ Y | **Z**. Si la relation est vraie, alors retirer l'arc correspondant;
4. procéder ainsi de suite avec des ensemble de taille 3, 4, etc.

Astuce: utilisez la fonction `combn(s, m)` pour obtenir toutes les combinaisons de taille `m` d'un ensemble `s` (attention le retour est sous forme de matrice). Pour obtenir le voisinage d'un noeud `x` dans un graphe `g`, utilisez `g$nodes[[x]]$nbr`. Pour supprimer un arc, utilisez `g = drop.edge(g, from = x, to = x)`. Pour une exécution rapide, préférez un seuil de tolérance faible (alpha=0.01).

## Orienter les arcs

Afin d'orienter les arcs, modifiez tout d'abord votre code de l'étape précédente afin de conserver chaque ensemble **Z**<sub>X,Y</sub> qui a permis de retirer un arc X-Y. Enfin, orientez le graphe en respectant la règle suivante:

1. orienter X -> W <- Y s'il n'y a pas d'arc entre X et Y et si W n'est pas inclus dans **Z**<sub>X,Y</sub>;
2. pour chaque arc non-orienté restant, décider une orientation arbitraire sans ajouter de nouvelle *v*-structure au graphe.

Astuce: utilisez `w %in% z` pour déterminer si un élément `w` est contenu dans un ensemble `z`. Pour placer un arc orienté, utilisez `g = set.arc(g, from = x, to = w)`.
