# Agsous-Monte-Carlo-program

# README - Utilisation du programme

## Description

Ce projet est un programme Python pour la modélisation de repliement protéiques utilisant un algorithme REMC (Replica Exchange Monte Carlo) et d'autres mouvements explicité dans l'article [A replica exchange Monte Carlo algorithm for protein folding in the HP model](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2071922/pdf/1471-2105-8-342.pdf). Il vous permet de générer une conformations optimal d'une protéineen 2D qui va minimiser l'energie globale de la protéine.

## Prérequis

- Python==3.10
- Conda==23.7.4

## Installation

Clonez ce référentiel sur votre machine :

   ```bash
   git clone git@github.com:dzSalim/Agsous-Monte-Carlo-program.git
```

Accéder au fichier du projet :

  ```bash
  cd Agsous-Monte-Carlo-program/
```

Créez l'environnement du projet :

   ```bash
   conda env create -f environment.yml
```

Activer l'environnement :

  ```bash
 conda activate Monte-Carlo-Project
```

# Lancer le programme

## Commande principale

```bash
python3 AGSOUS_Monte_Carlo_project.py [-h] (positional argument) PROTEIN.FASTA -n (optional argument) ITERATION
```

## Exemple d'utilisation 

```bash
python3 AGSOUS_Monte_Carlo_project.py P02776.fasta -n 1000
```

   
