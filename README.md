# pancreatic-cancer-subtypes-survival-linear-predictor
Show that there is no significant relationship between cancer subtypes and our linear predictor of survival



## DATA
4 sous types : 1squamous 2immunogenic 3progenitor 4ADEX
Récupérer en interrogeant la base GDC du TCGA (type = paad)
Fusion avec le prédicteur linéaire des modèles d’intérêt (8;9;10)

### Objectif : vérifier que les sous types de cancer sont indépendant du prédicteur linéaire continue (Anova) et dichotomiser sur Maxstat (Fisher) 

