Feature counts for each cell are divided by the total counts for that cell and multiplied by the `scale.factor`. This is then natural-log transformed using `log1p`
- Implementiert als Standard Normalize in Seurat.


## Erklärt anhand von Seurat-Implemtierung


data <- t(t(data) / Matrix::colSums(data)) * scale.factor
data <- log1p(data)
return(data)
scale.factor ist default 10000
Beispiel:
mat <- matrix(1:12, nrow = 3, ncol = 4)
print("Matrix:")
mat

[,1] [,2] [,3] [,4]
[1,] 1 4 7 10
[2,] 2 5 8 11
[3,] 3 6 9 12

gehen einfach von Ganzen Zahlen aus. Wir haben 3 gene und 4 Zellen.

t(t(mat) / Matrix::colSums(mat))

[,1][,2][,3][,4]
[1,] 0.1666667 0.2666667 0.2916667 0.3030303
[2,] 0.3333333 0.3333333 0.3333333 0.3333333
[3,] 0.5000000 0.4000000 0.3750000 0.3636364

t() steht für Transpose. **Zeile und Reihe werden vertauscht**. ColSums summiert alle Werte innerhalb einer
Spalte (z.b alle Werte in Spalte 1).Jeder Wert in der Reihe (vorher die Spalte) wird mit der **Summe der jeweiligen Spalte geteilt**. Es werden also eigentlich einzelne Zellenwerte (Spalten) von der Summe der Zellenwerte abgezogen (Summe der jeweiligen Spalte).

Data ← t(t(mat) / Matrix::colSums(mat)) * 10000

[,1] [,2] [,3] [,4]
[1,] 1666.667 2666.667 2916.667 3030.303
[2,] 3333.333 3333.333 3333.333 3333.333
[3,] 5000.000 4000.000 3750.000 3636.364

Als Default wird mit dem Wert 10.000 multipliziert.
log1p(data)

[,1] [,2] [,3] [,4]
[1,] 7.419181 7.888959 7.978539 8.016748
[2,] 8.112028 8.112028 8.112028 8.112028
[3,] 8.517393 8.294300 8.229778 8.199014

log1p() berechnet log(1+x) für alle Werte x in der Tabelle. Log() ist dabei der natürliche Log. Wird genutzt, da
Expressionswerte sehr klein sind anders als die Zahlen hier. Die +1 gewichtet sie.

TLDR Fassung:
Divide each cell by the total number of molecules measured in the cell
Multiply that number by a scaling factor (i.e. 10000)
Add 1, and take a natural log

## Vorteile von LogNormalize

Hi,
Usually we recommend using log-normalization because it has the advantages as followed:
1.Variance stabilization: scRNA-seq data is inherently noisy, and the counts can have high variability
due to technical factors and biological heterogeneity. Log-normalization helps stabilize the variance
across different expression levels. It reduces the influence of highly expressed genes and amplifies
the signal from lowly expressed genes, allowing for better identification of differentially expressed
genes and improved downstream analysis.
2.Interpretability: Log-normalization makes the data more interpretable by transforming the count
data into a log scale, which is more in line with how researchers usually think about fold-changes in
gene expression. This is especially useful when comparing gene expression across different cell
types or conditions, as fold-changes are more meaningful than absolute differences in counts.
3.Handling zero counts: A significant proportion of scRNA-seq data may contain zero counts for
some genes in certain cells due to dropout events (i.e., failure to detect lowly expressed genes). Log-
normalization typically adds a small pseudocount to the raw counts before taking the log
transformation, which helps handle zero counts and reduces the impact of dropout events on
downstream analysis.
Using "RC"only divide the count with a cell's total count, but without the log1p step as in log-normalization. It
is definitely okay to use RC, but you should do it based on the fact that it is appropriate for your data. You
should Not do it because it gives you more DEGs.