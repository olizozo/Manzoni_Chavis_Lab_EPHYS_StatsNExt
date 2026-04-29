# 🧬 Master Stats Suite : Analyse de Population
### *Chavis Lab & Manzoni Lab | Analyse de la Plasticité Synaptique*

Cette application est la station de travail finale du pipeline d'analyse. Elle permet de regrouper les données de plusieurs cellules/animaux, de visualiser les cinétiques de population et de réaliser des tests statistiques robustes (Modèles Mixtes et ANOVA) pour valider vos découvertes en plasticité synaptique (LTP/LTD).

👉 **Accès en ligne : [https://chavislab-plasticity.streamlit.app/](https://chavislab-plasticity.streamlit.app/)**

---

## 📋 Sommaire
1. [Fonctionnement de l'Interface](#1-fonctionnement-de-linterface)
2. [Détail des Analyses](#2-détail-des-analyses)
3. [Pourquoi ces Statistiques ?](#3-pourquoi-ces-statistiques)
4. [Limites & Bonnes Pratiques](#4-limites--bonnes-pratiques)
5. [Installation](#5-installation)

---

## 1. Fonctionnement de l'Interface

### La Barre Latérale (Sidebar)
* **Gestion des Groupes** : Vous pouvez ajouter autant de groupes que nécessaire (ex: WT, KO, KO+Drug).
* **Chargement CSV** : L'app accepte les fichiers CSV exportés par le script de biophysique.
* **Identifiant de Portée (Litter ID)** : **Crucial.** Vous devez renseigner quelle souris provient de quelle portée pour que le modèle mixte puisse corriger les biais génétiques familiaux.

### Onglet 1 : Visualisation & Export
* **Graphiques de Population** : Affiche la moyenne pondérée avec, au choix, l'erreur standard (SEM) ou l'intervalle de confiance à 95% (CI 95%).
* **Export Prism-Ready** : Génère un fichier au format "Wide" (large). Chaque colonne est déjà nommée (Mean, SEM, CI95) pour être collée directement dans GraphPad Prism.

---

## 2. Détail des Analyses

### Onglet 2 : Nested GLMM (Modèle Mixte)
Le Modèle Linéaire Mixte (LMM) est l'outil le plus puissant pour analyser la **maintenance** de la plasticité.
* **Fonction** : Il analyse la pente de la courbe (cinétique) après l'induction.
* **Interprétation** : Si l'interaction *Temps × Groupe* est significative ($p < 0.05$), cela prouve que la stabilité de la LTP/LTD est différente entre vos groupes.


### Onglet 3 : Stats Simples (Analyse par Époques)
Cette section découpe l'expérience en 3 blocs temporels (Baseline, Early, Late) à l'aide de curseurs.
1. **Test de Plasticité (vs 100%)** : Vérifie pour chaque groupe si la plasticité a été induite avec succès.
2. **ANOVA Inter-Groupes** : Vérifie s'il existe une différence globale entre les groupes sur une fenêtre précise.
3. **Taille d'Effet (Cohen's $d$)** : Quantifie la force du résultat biologique au-delà du simple $p < 0.05$.

---

## 3. Pourquoi ces Statistiques ?

### La Correction de Benjamini-Hochberg (FDR)
En testant la plasticité sur plusieurs groupes et plusieurs époques, on multiplie les chances de trouver un "faux positif" (Erreur de Type I). L'application applique automatiquement le **FDR** pour ajuster les p-values et garantir que vos découvertes ne sont pas dues au hasard.

### Le Post-Hoc de Tukey HSD
Lorsqu'une ANOVA à 3 groupes est significative, elle ne dit pas *qui* est différent de *qui*. Le test de Tukey compare toutes les paires possibles (ex: WT vs KO, WT vs Rescue, KO vs Rescue) tout en protégeant contre l'inflation de l'erreur globale.


---

## 4. Limites & Bonnes Pratiques

### Les Limites du Modèle (Convergence)
* **Le problème du N** : Le GLMM est complexe. Si vous n'avez qu'une seule portée ou seulement 2 animaux par groupe, le modèle risque de ne pas "converger" (échec du calcul). 
* **Action** : Dans ce cas, fiez-vous à l'onglet "Stats Simples" qui est plus robuste pour les petits échantillons.

### Bonnes Pratiques
* **Vérifiez la Baseline** : Assurez-vous que l'époque "Baseline" sélectionnée est bien stable à $100\%$.
* **Normalité** : Ces tests supposent une distribution relativement normale. Si vos données contiennent des "Outliers" extrêmes (ex: une cellule à $400\%$ alors que les autres sont à $120\%$), le modèle peut être faussé.

---

## 5. Installation

Pour les chercheurs souhaitant une version locale :

```bash
# 1. Installer les dépendances
pip install streamlit pandas numpy matplotlib statsmodels scipy

# 2. Lancer l'app
streamlit run app_nested_glmm.py
```

**Fichier `requirements.txt` requis pour le déploiement Cloud :**
```text
streamlit
pandas
numpy
matplotlib
statsmodels
scipy
```

---
*Développé pour le Chavis Lab | 2026; contact: olivier.manzoni@inserm.fr*
*Rigueur Statistique, Transparence Biologique.*
