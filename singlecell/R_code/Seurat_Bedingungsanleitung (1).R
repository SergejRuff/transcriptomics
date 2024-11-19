rm(list=ls())
# Seurat Anleitung

# Packages ----

library ( Seurat )
library ( patchwork )
library (ggplot2)
library (dplyr)
library(stringr)

# Import der Daten ----

list.files("GSE111429_RAW") # Die Dateien müssen folgende Namen tragen: "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"  

SC_RAW <- Read10X ( data.dir = "GSE111429_RAW" )

# Nutze Daten von Paper: https://www.nature.com/articles/s41467-020-14296-y#data-availability
# Vorteil: Sie liefern die originalen R-Skripte und mehr Daten für Reproduktion!



# erschaffe Seurat-Objekt ----

SeuratObj <- CreateSeuratObject(counts = SC_RAW, project = "Seurat_test_1", min.cells = 0, min.features = 0)

# min.cells = erste Filterung. Behalte nur Features, die mindestens in einer bestimmten Anzahl an Zellen expremiert wurden.
# min.Features = Behalte Zellen, wo eine mindest Anzahl an Features exprimiert wurde.
# Project = gebe jeder Zelle im Objekt einen Namen, welches die Zellen auf das Seurat_objekt zurueckfuehrt.
# Project Namen sind wichtig beim Mergen von Objekten, damit Zellen Patienten oder geweben zuweisen kann.


# Untersuche Elemente des Seurat Objektes ----

SeuratObj 

colnames ( x = SeuratObj ) # zeige alle Zellnamen
Cells ( SeuratObj ) # alternative
ncol ( SeuratObj )

rownames ( x = SeuratObj ) # zeige Featurenamen (gennamen)
Features(SeuratObj)
nrow ( SeuratObj ) #Anzahl an Features

# assay enthält die Rohdaten
colnames(SeuratObj@assays$RNA$counts) # Zellen in den Spalten
rownames(SeuratObj@assays$RNA$counts) # Gene in den Reihen

Idents ( SeuratObj ) # wichtig, wenn man Projekte vermischt.

SeuratObj_Metadaten <- SeuratObj@meta.data 
SeuratObj_Metadaten # zeige Metadaten der einzelnen ZELLEN. Wie viele Counts und NFeatures sie jeweils haben.
SeuratObj[[]] # Alternative Moglichkeit Metadaten aufzurufen

summary ( SeuratObj_Metadaten$nCount_RNA )
summary ( SeuratObj_Metadaten$nFeature_RNA )

str ( SeuratObj )
View ( SeuratObj ) # besser als str()

# assays. Alle Assays, die ich abgespeichert habe - bsp Rawdata und processed data
# meta.data: Alle Metadaten zu den Zellen. 
# active.assay: Welches Assay wird als default benutzt, wenn nichts anderes angegeben wird. DADRAUF ACHTEN!
#reductions: Welche Dimensionsreduktionen wurden abgespeichert
# graphs: Welche Plots (UMAP) wurden abgespeichert.

# https://satijalab.org/seurat/articles/essential_commands.html Seite mit Vignette fur mehr.

# Quality Control ----

# basierend auf: https://bioconductor.org/books/3.14/OSCA.basic/quality-control.html

# Motivation: Low quality libraries und die Probleme der Probenisolation beheben.

# Zellen mit hohem Mitochondrialen RNA-Anteil entfernen, da diese Zellen auf Zellschaden hinweisen.

sample_names <- grep("^mt", names(SeuratObj@assays$RNA$counts[,1]), value = TRUE)
sample_names
SeuratObj [[ "percent.mt" ]] <- PercentageFeatureSet ( SeuratObj, pattern = "^mt-" )
SeuratObj [[]] # haben eine neue Spalte in Meta.data hinzugefugt mit Anteil an Mitochondrialer RNA.
SeuratObj [[ "percent.mt" ]]
# Kontrolle Prozentzahl per Hand ausrechnen.
sum(SeuratObj@assays$RNA$counts[sample_names, 1])/sum(SeuratObj@assays$RNA$counts[, 1])*100

# Wichtig ! Habe eine eigene, eindeutige Endung fur mt-RNAs!

VlnPlot ( SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #VLNPlot

# im Vlnplot sehen wir wie hoch der Anteil an Zellen mit hohem mtRNA-Gehalt ist.
# Wir sehen wieviele Features die Zellen haben
# und wie viele Molekuele insgesamt.
# Zuviele oder zu wenige Molekuele und Features sprechen fur eine schlechte Qualitat und muessen rausgefiltert werden.
# Zu wenig: Ich habe Zellschaden und Dropouts, wodurch RNA-Anteil sinkt
# Zu viel: Mehrere Zellen hingen zusammen (nicht vernunftig in Suspension gebracht) und wurden zusammen sequenziert.

plot1 <- FeatureScatter(SeuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")+ geom_smooth(method="lm")
plot1 # Wenige nCounts aber hohe percent.mt -> Schlecht.
plot2 <- FeatureScatter ( SeuratObj, feature1 ="nCount_RNA", feature2 = "nFeature_RNA")+ geom_smooth(method="lm")
#plot1+plot2
plot2
# mit Featurescatter plotte ich zwei variablen gegeneinander.
# Gute Qualitat = sollte der Linie folgen.

SeuratObj <- subset ( SeuratObj, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & percent.mt < 5 )
# Behalte nur Zellen mit mehr als 1500 und weniger als 5000 Features und weniger als 5 % mt Anteil.

# plot again
VlnPlot ( SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #VLNPlot
plot3 <- FeatureScatter ( SeuratObj, feature1 ="nCount_RNA", feature2 = "nFeature_RNA")+ geom_smooth(method="lm")
plot2 +plot3
# Daten normalisieren ----

#  Um Genexpression über mehrere Zellen zu vergleichen, muss man Daten normalisieren.

SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObj[["RNA"]]@layers$data # aufrufen der normalisierten Werte. Daten in Seurat Objekt gespeichert.

# Feature Selection ----

# Konzentrieren uns nur auf Variable Gene,da diese Zelltypisch sind. Keine Housekeeping gene.

SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
# identifiziere hochvariable Features.

top10 <- head(VariableFeatures(SeuratObj), 10) # zeige top 10 variable features
top10

plot1 <- VariableFeaturePlot(SeuratObj)
# Add labels of top10 genes
plot2 <- LabelPoints(plot = plot1, points = top10, repel=TRUE, xnudge = 0, ynudge=0)
plot2

#Daten skalieren ----

# Entzferne Varianz/ Noice, die auf technische Fehler basieren.Batch_Effekte.

# Wichtiger Schritt fur die PCA.

all.genes <- rownames(SeuratObj)
SeuratObj <- ScaleData(SeuratObj, features = all.genes)

# Skalierung ist nowendig, damit hoch exprimierte Gene nicht die weiteren Schritte dominieren.
# Mean expression wird auf 0 gesetzt
# Varianz zwischen allen Zellen auf 1 gesetzt.

SeuratObj@assays[["RNA"]]@layers[["scale.data"]] # Werte der Scaladata aufrufen

# PCA ----

# Veruschen so Quellen der Variabilität zu finden, die Daten am besten erklären.

# braucht skalierte Daten + variable Features.
SeuratObj <- RunPCA( SeuratObj , features = VariableFeatures(object = SeuratObj)) 

colnames(SeuratObj@reductions[["pca"]]) # pca dim
rownames(SeuratObj@reductions[["pca"]]) # cell

print(SeuratObj[["pca"]], dims = 1:5, nfeatures = 5) # hoch vs runterreguliert.

#DimPlot(SeuratObj, reduction ="pca")

#SeuratObj <- JackStraw(SeuratObj, num.replicate = 100)
#SeuratObj <- ScoreJackStraw(SeuratObj, dims = 1:20)

#JackStrawPlot(SeuratObj, dims = 1:15)
# Jackstaw - Plottet p values für alle PCA Dim.

# behalte nur PC, die relevant sind. Müssen varianz erklären.
ElbowPlot(SeuratObj) # subjektiv

# Zellen clustern ----

SeuratObj <- FindNeighbors(SeuratObj, dims = 1:15) # finde Nachbarn für Zellen.
SeuratObj <- FindClusters(SeuratObj, resolution = 0.5)
# Je höher die Resolution, desto mehr Cluster werden gefunden.

head(Idents(SeuratObj), 5) # zeige an, in welche Cluster die ersten 5 Zellen gehoeren.

# UMAP ----

SeuratObj<- RunUMAP(SeuratObj, dims = 1:10)
DimPlot(SeuratObj, reduction = "umap",label=TRUE)

# Cluster finden (welche Gene definieren Cluster) ----

# Finde DEG pro Cluster. Bis jetzt nur Variable Gene für gesamten Datensatz untersucht.

cluster2.markers <- FindMarkers(SeuratObj, ident.1 = 2, min.pct = 0.25) # einzelne Cluster untersuchen
head(cluster2.markers, n = 5)

# Alle Cluster untersuchen 

# ?FindConservedMarkers -> besser, wenn man zwei gruppen vergleicht.

pbmc.markers <- FindAllMarkers(SeuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(SeuratObj, features = c("Hcar2", "Cxcl14", "Serpina3g")) #0,0,4

DimPlot(SeuratObj, reduction = "umap",label=TRUE)

# rename cluster

SeuratObj <- RenameIdents(SeuratObj,"3"="test_zellen")
Idents(SeuratObj)

# Daten integrieren/verbinden ----

#merge(erstes_Seurat_Object, y= (zweites_seuratObj , drittes_SeuratObj, viertes_SeuratObj),
#      add.cell.ids = ls()[]) # add.cell.ids ist sollte nicht notwendig sein, wenn du gute Projectnamen hast.

# Danach wurde man die weiteren Schritte ab Quality Control durchfuhren

# Weitere Quellen ----

# https://bioconductor.org/books/3.14/OSCA.basic/quality-control.html
# https://satijalab.org/seurat/archive/v3.0/integration.html
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://www.youtube.com/watch?v=5HBzgsz8qyk&t=461s
#https://www.youtube.com/watch?v=HrbeaEJqKcY
#https://satijalab.org/seurat/archive/v3.0/pancreas_integration_label_transfer.html tutorial fur integration
# https://www.youtube.com/watch?v=PSuiBaUqh0Y&t=337s Video zum Thema Integration