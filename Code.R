annotation<-read.table("GO_screen_annot.txt", header=T, sep="\t")
head(annotation)
library(Biobase)
micro_eset<-readRDS("micro_eset.RDS")

fdata_micro_eset <- fData(micro_eset)
head(fdata_micro_eset)

expr_micro_eset <- exprs(micro_eset)
expr_micro_eset

pdata_micro_eset <- pData(micro_eset)
pdata_micro_eset

pdata_micro_eset$ApoE

library(dplyr)
pdata_wt_micro_eset <- pdata_micro_eset %>%  filter (Genotype != "CC")

pdata_wt_micro_eset

expr_wt <- expr_micro_eset[, colnames(expr_micro_eset) %in% rownames(pdata_wt_micro_eset)]
head(expr_wt)
expr_wt0<-merge(expr_wt,annotation,by.x="row.names",by.y="transcript_cluster_id")
head(expr_wt0)
rownames(expr_wt0)<-make.names(expr_wt0$SYMBOL,unique=T)
expr_wt0 <- expr_wt0[, !(colnames(expr_wt0) %in% c("gene", "SYMBOL", "Row.names"))]
head(expr_wt0)

library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2.5,
cex.axis = 1.5, cex.main = 2)

clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

traitData <- pdata_wt_micro_eset
head(traitData)
traitData$ApoE.Treatment<-paste(traitData$ApoE, traitData$Treatment, sep=".")
table(traitData$ApoE.Treatment)

traitData$Genotype.Treatment <- NULL

traitData$Untreated <- ifelse(traitData$Treatment == "Untreated", 1, 0)  
traitData$Myelin <- ifelse(traitData$Treatment == "Myelin", 1, 0)  
traitData$LPS <- ifelse(traitData$Treatment == "LPS", 1, 0)

traitData$APOE22<-ifelse(traitData$ApoE == "APOE22", 1, 0)
traitData$APOE33<-ifelse(traitData$ApoE == "APOE33", 1, 0)
traitData$APOE44<-ifelse(traitData$ApoE == "APOE44", 1, 0)
traitData$KO<-ifelse(traitData$ApoE == "KO", 1, 0)

traitData$APOE22.LPS <- ifelse(traitData$ApoE.Treatment == "APOE22.LPS", 1, 0)
traitData$APOE33.LPS <- ifelse(traitData$ApoE.Treatment == "APOE33.LPS", 1, 0)
traitData$APOE44.LPS <- ifelse(traitData$ApoE.Treatment == "APOE44.LPS", 1, 0)
traitData$KO.LPS <- ifelse(traitData$ApoE.Treatment == "KO.LPS", 1, 0)

traitData$APOE22.Myelin <- ifelse(traitData$ApoE.Treatment == "APOE22.Myelin", 1, 0)
traitData$APOE33.Myelin <- ifelse(traitData$ApoE.Treatment == "APOE33.Myelin", 1, 0)
traitData$APOE44.Myelin <- ifelse(traitData$ApoE.Treatment == "APOE44.Myelin", 1, 0)
traitData$KO.Myelin <- ifelse(traitData$ApoE.Treatment == "KO.Myelin", 1, 0)

traitData$APOE22.Untreated <- ifelse(traitData$ApoE.Treatment == "APOE22.Untreated", 1, 0)
traitData$APOE33.Untreated <- ifelse(traitData$ApoE.Treatment == "APOE33.Untreated", 1, 0)
traitData$APOE44.Untreated <- ifelse(traitData$ApoE.Treatment == "APOE44.Untreated", 1, 0)
traitData$KO.Untreated <- ifelse(traitData$ApoE.Treatment == "KO.Untreated", 1, 0)

table(traitData$LPS)

datTraits <- traitData[rownames(datExpr), , drop = FALSE]

datTraits <- datTraits[, c("Untreated", "Myelin", "LPS", "KO",
                           "APOE22", "APOE33", "APOE44",
                           "APOE22.LPS", "APOE33.LPS", "APOE44.LPS","KO.LPS",
                           "APOE22.Myelin", "APOE33.Myelin", "APOE44.Myelin","KO.Myelin",
                           "APOE22.Untreated", "APOE33.Untreated", "APOE44.Untreated", "KO.Untreated")]

print(head(datTraits))


traitData_filter = traitData[15:33]
table(traitData$Treatment)  
table(traitData$LPS)


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData_filter, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(traitData_filter),
main = "Sample dendrogram and trait heatmap")

table(traitData$Treatment)  # Verifica si "LPS" está en el tratamiento
table(traitData$LPS)        # Verifica si hay valores 0 y 1 en 'traitData$LPS'

powers = c(c(1:10), seq(from = 12, to=16, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr0, power = 12,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "microglia",
verbose = 3)


mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],  
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

table(mergedColors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "Microglia-02-iPSC.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData_filter, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(traitData_filter),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 1,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

head(moduleTraitCor)
head(moduleTraitPvalue)
# head(textMatrix)

head(cbind(data.frame(moduleTraitCor[,1]),data.frame(moduleTraitPvalue[,1])))

dfcor<-cbind(data.frame(moduleTraitCor[,1]),data.frame(moduleTraitPvalue[,1]))
colnames(dfcor)<-c("R_KO","P_KO")
head(dfcor)
write.table(dfcor,"module.correlations.txt",sep="\t")

dfcor[order(dfcor$P_KO,decreasing=F),]

Untreated = as.data.frame(traitData$Untreated)
Myelin = as.data.frame(traitData$Myelin)
LPS = as.data.frame(traitData$LPS)
KO = as.data.frame(traitData$KO)
APOE22 = as.data.frame(traitData$APOE22)
APOE33 = as.data.frame(traitData$APOE33)
APOE44 = as.data.frame(traitData$APOE44)
APOE22.LPS = as.data.frame(traitData$APOE22.LPS)
APOE33.LPS = as.data.frame(traitData$APOE33.LPS)
APOE44.LPS = as.data.frame(traitData$APOE44.LPS)
KO.LPS = as.data.frame(traitData$KO.LPS)
APOE22.Myelin = as.data.frame(traitData$APOE22.Myelin)
APOE33.Myelin = as.data.frame(traitData$APOE33.Myelin)
APOE44.Myelin = as.data.frame(traitData$APOE44.Myelin)
KO.Myelin = as.data.frame(traitData$KO.Myelin)
APOE22.Untreated = as.data.frame(traitData$APOE22.Untreated)
APOE33.Untreated = as.data.frame(traitData$APOE33.Untreated)
APOE44.Untreated = as.data.frame(traitData$APOE44.Untreated)
KO.Untreated = as.data.frame(traitData$KO.Untreated)

# Calcular la correlación de cada gen con los módulos
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))

# Poner nombres a las columnas
names(geneModuleMembership) = paste("MM", names(MEs), sep="")

# Definir las variables a partir de traitData
variables = list(Untreated, Myelin, LPS, KO, APOE22, APOE33, APOE44, 
                 APOE22.LPS, APOE33.LPS, APOE44.LPS, KO.LPS, 
                 APOE22.Myelin, APOE33.Myelin, APOE44.Myelin, KO.Myelin, 
                 APOE22.Untreated, APOE33.Untreated, APOE44.Untreated, KO.Untreated)

var_names = c("Untreated", "Myelin", "LPS", "KO", "APOE22", "APOE33", "APOE44", 
              "APOE22_LPS", "APOE33_LPS", "APOE44_LPS", "KO_LPS", 
              "APOE22_Myelin", "APOE33_Myelin", "APOE44_Myelin", "KO_Myelin", 
              "APOE22_Untreated", "APOE33_Untreated", "APOE44_Untreated", "KO_Untreated")

# Definir los nombres de los módulos (modNames) y los colores de los módulos (moduleColors)
modNames = substring(names(MEs), 3)  # Nombres de módulos
moduleColors = as.vector(moduleColors)  # Colores de los módulos
modules = c("blue", "brown", "turquoise", "grey","red","black","yellow","green")  # Módulos específicos a evaluar

# Función para calcular las correlaciones y realizar gráficos de dispersión para cada variable
plot_gene_module_correlation <- function(var_name, var_data, geneModuleMembership, geneTraitSignificance, moduleColors) {
  
  # Calcular la correlación para esta variable
  geneTraitSignificance = as.data.frame(cor(datExpr0, var_data, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  # Nombre de las columnas
  names(geneTraitSignificance) = paste("GS.", var_name, sep="")
  names(GSPvalue) = paste("p.GS.", var_name, sep="")
  
  # Graficar los resultados para cada módulo
  for (module in modules) {
    column = match(module, modNames)  # Seleccionar la columna correspondiente al módulo
    moduleGenes = moduleColors == module  # Filtrar genes del módulo actual
    
    # Crear gráfico
    par(mfrow = c(1, 1))  # Configurar el gráfico
    verboseScatterplot(
      abs(geneModuleMembership[moduleGenes, column]),  # Membresía del módulo
      abs(geneTraitSignificance[moduleGenes, 1]),  # Significancia del gen para la variable
      xlab = paste("Module Membership in", module, "module"),
      ylab = paste("Gene significance for", var_name),
      main = paste("Module membership vs. gene significance for", var_name, "(", module, "module)"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module  # Colorear según el módulo
    )
  }
}

# Ejecutar para cada variable
for (i in 1:length(variables)) {
  var_name = var_names[i]
  var_data = variables[[i]]  # Tomar los datos específicos de la variable
  
  # Llamar a la función para calcular la correlación y generar los gráficos
  plot_gene_module_correlation(var_name, var_data, geneModuleMembership, geneTraitSignificance, moduleColors)
}


# Crear la carpeta "plots" si no existe
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Función para calcular las correlaciones y generar gráficos de dispersión guardados
plot_gene_module_correlation <- function(var_name, var_data, geneModuleMembership, geneTraitSignificance, moduleColors) {
  
  # Calcular la correlación para esta variable
  geneTraitSignificance = as.data.frame(cor(datExpr0, var_data, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  # Nombre de las columnas
  names(geneTraitSignificance) = paste("GS.", var_name, sep="")
  names(GSPvalue) = paste("p.GS.", var_name, sep="")
  
  # Guardar los gráficos para cada módulo
  for (module in modules) {
    column = match(module, modNames)  # Seleccionar la columna correspondiente al módulo
    moduleGenes = moduleColors == module  # Filtrar genes del módulo actual
    
    # Crear nombre del archivo de salida para el gráfico
    output_file = paste0("plots/", var_name, "_", module, "_scatterplot.png")
    
    # Guardar el gráfico como imagen
    png(output_file, width = 800, height = 600)  # Guardar como PNG, ajusta las dimensiones
    verboseScatterplot(
      abs(geneModuleMembership[moduleGenes, column]),  # Membresía del módulo
      abs(geneTraitSignificance[moduleGenes, 1]),  # Significancia del gen para la variable
      xlab = paste("Module Membership in", module, "module"),
      ylab = paste("Gene significance for", var_name),
      main = paste("Module membership vs. gene significance for", var_name, "(", module, "module)"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module  # Colorear según el módulo
    )
    dev.off()  # Cerrar el dispositivo gráfico después de guardar
  }
}

# Ejecutar para cada variable
for (i in 1:length(variables)) {
  var_name = var_names[i]
  var_data = variables[[i]]  # Tomar los datos específicos de la variable
  
  # Llamar a la función para calcular la correlación y generar los gráficos
  plot_gene_module_correlation(var_name, var_data, geneModuleMembership, geneTraitSignificance, moduleColors)
}


# Hubs

ADJ1=abs(cor(datExpr0,use="p"))^15
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1[order(-Alldegrees1$kTotal) ,], n=20)
b<-chooseTopHubInEachModule(
   datExpr0, 
   moduleColors, 
   omitColors = NA, 
   power = 10, 
   type = "unsigned")
b
b<-chooseTopHubInEachModule(
   datExpr0, 
   moduleColors, 
   omitColors = NA, 
   power = 10, 
   type = "signed")
b

# Enrichment Analysis

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

# Crear el data frame inicial con la información de KO
geneTraitSignificance_KO = as.data.frame(cor(datExpr0, traitData$KO, use = "p"))
GSPvalue_KO = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_KO), nSamples))

# Asignar nombres a las columnas
names(geneTraitSignificance_KO) = "GS.KO"
names(GSPvalue_KO) = "p.GS.KO"

# Crear geneInfo0 con los datos iniciales
geneInfo0 = data.frame(
  moduleColor = moduleColors,
  geneTraitSignificance_KO,
  GSPvalue_KO
)

# Añadir significancia para cada una de las otras variables
for (i in 1:length(variables)) {
  var_name = var_names[i]
  var_data = variables[[i]]

  # Calcular correlación y p-valor
  geneTraitSignificance = as.data.frame(cor(datExpr0, var_data, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

  # Asignar nombres adecuados
  names(geneTraitSignificance) = paste("GS.", var_name, sep = "")
  names(GSPvalue) = paste("p.GS.", var_name, sep = "")

  # Añadir al data frame geneInfo0
  geneInfo0 = cbind(geneInfo0, geneTraitSignificance, GSPvalue)
}

# Ordenar los módulos según su correlación con KO
modOrder = order(-abs(cor(MEs, traitData$KO, use = "p")))

# Añadir la información de membresía en módulos
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(
    geneInfo0, 
    geneModuleMembership[, modOrder[mod]], 
    MMPvalue[, modOrder[mod]]
  )
  names(geneInfo0) = c(
    oldNames, 
    paste("MM.", modNames[modOrder[mod]], sep = ""), 
    paste("p.MM.", modNames[modOrder[mod]], sep = "")
  )
}

# Ordenar los genes primero por módulo, luego por significancia con KO
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.KO))
geneInfo = geneInfo0[geneOrder, ]

# Mostrar el resultado
head(geneInfo)



library(org.Hs.eg.db, lib = "/home/julian/R/library")
selected<-select(org.Hs.eg.db, as.character(rownames(geneInfo)), "ENTREZID", "SYMBOL")
head(selected)
selected2<-merge(selected,geneInfo,by.x="SYMBOL",by.y="row.names")
head(selected2)
selected2[which(selected2$SYMBOL=="APOE"),]

moduleColors<-selected2$moduleColor

GOenr = GOenrichmentAnalysis(moduleColors, selected2$ENTREZID, organism = "human", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.txt", sep = "\t", quote = TRUE, row.names = FALSE)
head(tab)
library(dplyr)
tab.table <-tab %>%
    group_by(module) %>%
    slice_max(n = 5, order_by = enrichmentP)
# tab.table[,c(1,5,6,13)]

tab.table[which(grepl("blue|brown|turquoise|grey|yellow|red|green|black",tab.table$module)),c(1,5,6,13)]

blue<-names(datExpr)[moduleColors=="blue"]
brown<-names(datExpr)[moduleColors=="brown"]
turquoise<-names(datExpr)[moduleColors=="turquoise"]
grey<-names(datExpr)[moduleColors=="grey"]
yellow<-names(datExpr)[moduleColors=="yellow"]
red<-names(datExpr)[moduleColors=="red"]
green<-names(datExpr)[moduleColors=="green"]
black<-names(datExpr)[moduleColors=="black"]


print("blue")
length(blue)
sort(blue)
print("brown")
length(brown)
sort(brown)
print("turquoise")
length(turquoise)
sort(turquoise)
print("grey")
length(grey)
sort(grey)

print("yellow")
length(yellow)
sort(yellow)

print("red")
length(red)
sort(red)

print("green")
length(green)
sort(green)

print("black")
length(black)
sort(black)

library(gprofiler2, lib = "/home/julian/R/library")

blue.Enrich<-gost(blue, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfblue<-as.data.frame(blue.Enrich$result)
dfblue <- data.frame(apply(dfblue,2,as.character))
head(dfblue[,c(3,10,11)])

dfblue.table <-dfblue %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfblue.table[,c(3,10,11)]

print("**blue**")
pblue <- gostplot(blue.Enrich, capped = FALSE, interactive = TRUE) 
pblue

write.table(dfblue,"blue.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2, lib = "/home/julian/R/library")

brown.Enrich<-gost(brown, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfbrown<-as.data.frame(brown.Enrich$result)
dfbrown <- data.frame(apply(dfbrown,2,as.character))
head(dfbrown[,c(3,10,11)])

dfbrown.table <-dfbrown %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfbrown.table[,c(3,10,11)]

print("**brown**")
pbrown <- gostplot(brown.Enrich, capped = FALSE, interactive = TRUE) 
pbrown

write.table(dfbrown,"brown.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2, lib = "/home/julian/R/library")

turquoise.Enrich<-gost(turquoise, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfturquoise <-as.data.frame(turquoise.Enrich$result)
dfturquoise <- data.frame(apply(dfturquoise,2,as.character))
head(dfturquoise[,c(3,10,11)])

dfturquoise.table <-dfturquoise %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfturquoise.table[,c(3,10,11)]

print("**turquoise**")
pturquoise<- gostplot(turquoise.Enrich, capped = FALSE, interactive = TRUE) 
pturquoise

write.table(dfturquoise,"turquoise.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2, lib = "/home/julian/R/library")

grey.Enrich<-gost(grey, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfgrey <-as.data.frame(grey.Enrich$result)
dfgrey<- data.frame(apply(dfgrey,2,as.character))
head(dfgrey[,c(3,10,11)])

dfgrey.table <-dfgrey %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfgrey.table[,c(3,10,11)]

print("**grey**")
pgrey <- gostplot(grey.Enrich, capped = FALSE, interactive = TRUE) 
pgrey

write.table(dfgrey,"grey.Enrichments.txt", sep="\t", row.names=F)


library(gprofiler2, lib = "/home/julian/R/library")

yellow.Enrich<-gost(yellow, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfyellow<-as.data.frame(yellow.Enrich$result)
dfyellow <- data.frame(apply(dfyellow,2,as.character))
head(dfyellow[,c(3,10,11)])

dfyellow.table <-dfyellow %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfyellow.table[,c(3,10,11)]

print("**yellow**")
pyellow <- gostplot(yellow.Enrich, capped = FALSE, interactive = TRUE) 
pyellow

write.table(dfyellow,"yellow.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2, lib = "/home/julian/R/library")

red.Enrich<-gost(red, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfred <-as.data.frame(red.Enrich$result)
dfred <- data.frame(apply(dfred,2,as.character))
head(dfred[,c(3,10,11)])

dfred.table <-dfred %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfred.table[,c(3,10,11)]

print("**red**")
pred <- gostplot(red.Enrich, capped = FALSE, interactive = TRUE) 
pred

write.table(dfred,"red.Enrichments.txt", sep="\t", row.names=F)

library(gprofiler2, lib = "/home/julian/R/library")

green.Enrich<-gost(green, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfgreen <-as.data.frame(green.Enrich$result)
dfgreen <- data.frame(apply(dfgreen,2,as.character))
head(dfgreen[,c(3,10,11)])

dfgreen.table <-dfgreen %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfgreen.table[,c(3,10,11)]

print("**green**")
pgreen <- gostplot(green.Enrich, capped = FALSE, interactive = TRUE) 
pgreen

write.table(dfgreen,"green.Enrichments.txt", sep="\t", row.names=F)


library(gprofiler2, lib = "/home/julian/R/library")

black.Enrich<-gost(black, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr",
domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)

dfblack<-as.data.frame(black.Enrich$result)
dfblack <- data.frame(apply(dfblack,2,as.character))
head(dfblack[,c(3,10,11)])

dfblack.table <-dfblack %>%
    group_by(source) %>%
    slice_max(n = 5, order_by = p_value)
dfblack.table[,c(3,10,11)]

print("**black**")
pblack <- gostplot(black.Enrich, capped = FALSE, interactive = TRUE) 
pblack

write.table(dfblack,"black.Enrichments.txt", sep="\t", row.names=F)


