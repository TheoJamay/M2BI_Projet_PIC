M2BI_Projet_PIC  
# Projet Protein Interactions Calculator (PIC)

Protein Interactions Calculator (PIC) est un serveur web qui reconnaît différents types d'interactions, telles que les ponts disulfure, les interactions hydrophobes, les interactions ioniques, les liaisons hydrogène, les interactions aromatique-aromatique, les interactions aromatique-sulfure et les interactions cation-pi au sein d'une même chaîne protéique ou entre chaîne protéines différente. Le fichier d'entrée doit être au format Protein data bank(.pdb).

## Requirement  

* Python 3.x (>= 3.6)  
* Python modules: NumPy, Math, Sys.  
* Reduce : https://github.com/rlabduke/reduce

## Running  

Le scrip s'exécute obligatoirement avec un fichier .pdb préalablement traité avec l'outil Reduce. Une précision importante est que le fichier traité par reduce doit être nommé de la manière suivante : myfileFH.pdb pour que le script puisse fonctionner.

```reduce -FLIP myfile.pdb > myfileFH.pdb```  
```python Intraprotein_interaction_PIC.py myfileFH.pdb```

## Output

Le script Intraprotein_interaction_PIC.py produit 8 fichiers de sortie au format .txt. Soit un fichier de sortie pour chaque type d'intéractions.
