
## Intro

Hier werden die verschiedenen Algos erklärt, die benutzt werden, um Gendaten für ScRNA-Seq zu normalisieren. Verschiedene Beispiele werden gegeben  und Sie werden erklärt anhand von konkreten Code-Beispielen (in Seurat z.B.)

## Wozu braucht man Normalisierung?


## Beispiele für Algos


- [[LogNormalize]]: Feature counts for each cell are divided by the total counts for that cell and multiplied by the `scale.factor`. This is then natural-log transformed using `log1p`