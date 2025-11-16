#############################################################################
#
# PRACTICA R
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################
###primero descargamos el paquete y donde pone URL le pegamos la direccion web de los datos, y donde pone data "read table" leemos los archivos de ese enlace sin necesidad de desacrar los archivos####
# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data)
head(data)
tail(data)
###dim: 12488 genes que tenemos en 6 columnas###
## head data: las 6 primeras filas de esta tabala y nos sirve para ver como estan organizados, ver como se llaman las columnas WT: wild type y KO: knockout tenemos 3 replicacas para cada condicion
##tail data: nos salen las 6 ultimas filas 
# Hacemos un primer histograma para explorar los datos
hist(data, col = "gray", main="GSE5583 - Histogram")
# abrir el data en histograma asi a lo bruto, el grey para el color y el main como para el titulo

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve? para que quede bonito y estos datos los guardamos en data2 y son los datos transformados
data2 = log2(data)
hist(data2, col = "gray", main="GSE5583 (log2) - Histogram")
# esto seria la transformacion de los datos, y para las fotos le ponemos mas cosas para ponerlo bonito 

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot? hay 1 por cada columna, azul los WT y naranja los KO, la linea negra es la mediana osea Q2, arriba el Q1 y abajo el Q3, los bigotes la disperción, los puntos en los extremos son los valores extremos 
boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)
#que es "las=2" es para el nombre de las columnas, con salen los nombres verticales y sin pues salen horizontal y no salen todos	

#usamos los datos transformados con log, donde pone col para los colores, se ponen en orden del mismo orden de la tabla, titulo y fin 
boxplot(data, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)
# si le cambiamos al data normal queda mas feo pq no sale 

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
# de los valores de expresión. ¿Es correcta la separación?
hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")
#que es, para ver si lo hemos hecho bien nos pone los datos en valores de expresion osea si tenemos 3 casos control deberian de ser similares y si tenemos 3 de enfermeos deberian de clasificarse las 3 juntitas y luego estas dos separadas, tenemos uno en el ejemplo y una rama de WT y otra de KO entonces que lo hemos hecho bien si una WT nos sale en la rama de Ko es que hay una que nos esta dando error 


#######################################
# Análisis de Expresión Diferencial 
#######################################
head(data)
#vamos a coger las WT y las KO y las vamos a comparar e ir viendo gen a gen

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado? matriz 
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)
#dentro del corchete para buscar en tabla antes de la coma son filas y detras de la coma las columnas y si no ponemos nada delante son todas las filas
head(wt)
head(ko)
# vemos los datos de cada columna, class para ver el tipo de objeto 

# Calcula las medias de las muestras para cada condición. Usa apply sacamos la media en wt y ko, calculamos con función apply y como funciona, el apply con el objeto que queremos calcular un 1 para todas las filas osea para cada gen o 2 para las columnas y mean que es la media sale arriba el gen y abajo la media

wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)

# ¿Cuál es la media más alta? 37460.5 
limit = max(wt.mean, ko.mean)
limit

# Ahora hacemos un scatter plot (gráfico de dispersión) es un garfico de dispercion en x wt y en y el ko, la lista de media de ko con el - fancy este el wt, que es xlab, la etiqueta del x xlim e ylim el limite
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Añadir una línea diagonal con abline
abline(0, 1, col = "red")
# la linea nos indicaria una correlacion perfecta en concreto la funcion abline es para añadir lineas 

# ¿Eres capaz de añadirle un grid? si solo darle command enter 
grid()
#abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
abline(1, 2, col = "pink")     # línea y = 2x + 1
abline(h = 2, col = "green")  # línea y = 2
abline(v = 3, col = "violet") # línea x = 3

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?
pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
# los puntos son los genes, el eje x pone diferencia de medias, en el caso del 3 el corte en 0 implica que no hay cambio despues indica cambios en la exprecion si hay sobrexprecion de ko serian negativos izq caso contrario osea exprecion wt ademas que estan reprimidos en ko, linea verde corte de pvalue como esta con el log lo pone en el numero 2 y los significativos nos salen por arriba de la verde y por la azul 

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)
#cuantos genes nos salen que pasan ambos genes, los genes son las filas y nos queda 426 que son diferencialmente expresados 

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")
# las pelotitas rojas son genes sobreexpresados en ko y los azules son reprimidos en ko o los sobreexpresados en wt y los negros son los que no pasan el filtro

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)
#es un mapa de calor nos permite combrobar la clasificacion de la muestra y que estan bien clasificado los wt juntos y los ko juntos, se tienen que hacer con los genes exprecionalmente diferenciados si esta sobreexpresado en ko es pq estan reprimidos en wt lo que se ve en el heat map con la diferencia de color
heatmap(filtered)
#


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")	
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")
#rojo sobreexprecion y azul reprimido

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE)
