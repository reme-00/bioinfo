# Como en la seccion de la practica del aula virtual no está el archivo 
# del guion de prácticas en fomrato .qmd, he creado un archivo nuevo sobre el 
# que escribiré los códigos. 

library(phangorn)
library(phangorn)
mtdna <- read.phyDat('../../data/mtDNA.fasta', format = 'fasta')
image(mtdna)
glance(mtdna)
as.StringSet(mtdna)

#EJERCICIO: ¿Puedes confirmar el número de secuencias? ¿Qué longitud tiene el 
#alineamiento? ¿Qué función usarías para extraer los nombres de las secuencias?

#numero de secuencias:
length(mtdna)

#longitud del alineamiento:
attr(mtdna, "nr")

#extraccion de los nombres de las secuencias: 
names(mtdna)

#calculo de distancias: 
nj(dist.ml(mtdna))

#otra opcion (de phangorn): 
dist.ml(mtdna)

distancias <- dist.ml(mtdna, model = 'F81')
mt_nj    <- NJ(distancias)
mt_upgma <- upgma(distancias)

#EJERCICIO: 
#¿De qué clase son los objetos mt_nj y mt_upgma?
    #Lo que hacemos con ello es representar distancias geneticas entre hominidos
    #ambos pertenecen a la clase de phylo
#Consulta la estructura de estos objetos (str(mt_nj)).
str(mt_nj)
    #el primer $indica que se conectan los nodos, el segundo es un vector con la 
    #longitud de cada rama, el tercero son los nombres de las 42 secuencias y el 
    # cuarto la cantidad de nodos internos. 
#Represéntalos gráficamente (plot(mt_nj)) y consulta la ayuda de la función 
 #plot.phylo() para saber qué opciones de representación gráfica tienes.
plot(mt_nj)
plot.phylo(mt_nj)
?plot.phylo
#para hacerlo en formato circular
plot(mt_nj, type = "fan", cex = 0.7) 
#Para qué sirven las funciones add.scale.bar() y axisPhylo()?
add.scale.bar(mt_nj)
axisPhylo(mt_nj)
    #sirve, la primera para añadir una barra de escala en el gráfico, y la 
    #segunda para añadir un eje graduado. 
#vamos a comprobar que el arbol no está enraizado: 
is.rooted(mt_nj)
#opciones para enraizarlo: (debemos especificar el "outgroup" o el nodo) 
mt_nj_rooted <- root(mt_nj, node = 75)
#para identificar los nodos internos: 
nodelabels(mt_nj)
#¿no entiendo por qué es el nodo 75 el ancestro entre el gorila y el chimpancé?
#¿como se ve eso?
#para asegurar que el enraizado no depende de la numeracion: 
# Debería dar el mismo resultado:
mt_nj_rooted <- root(mt_nj, outgroup = 'G.gorilla_D38114')
#arbol enraizado: 
plot(mt_nj_rooted)
#comando para escoger la posicion de la raiz con la que el arbol esté equilibrado
midpoint(mt_nj)
# Esto sobreescribe el objeto anterior:
mt_nj_rooted <- midpoint(mt_nj)
plot(mt_nj_rooted, align.tip.label = TRUE, main = 'NJ, enraizado en punto medio')
nodelabels(frame = 'none')
add.scale.bar()
#metodo de enraizado root-to-tip:
#-1º creamos un vector de fechas, pero primero debemos saber el orden en que 
#debemos especificar las edades: 
# 1. Asegúrate de que 'hojas' existe primero
hojas <- mt_nj$tip.label

# 2. Crea el vector (aquí es donde se define el nombre 'fechas')
fechas <- rep(0, length(hojas))

# 3. Asigna los nombres a ese vector
names(fechas) <- hojas
fechas <- c(0, 0, 0, 0, 0, 0, 0, 0, 0,
            -400000,
            -39000, -50000, -110000, -110000,
            -39820, -40096, -38515, -38310, -39000, -40000, -41210,
            -42430, -42540, -43230, -43780, -44290, -44770, -49000,
            -50000, -82752, -65000, -122287, -110450,
            0, 0, 0, 0, 0, 0, 0, 0, 0)
# Observa que hay que asignar el resultado de setNames() a "fechas"
fechas <- setNames(fechas, mt_nj$tip.label)
fechas
#posicion de la raiz que optimiza la regresion: 
mt_nj_rtt <- rtt(mt_nj, tip.dates = fechas, objective = 'rms')

#EJERCICIO: 
#¿Se corresponden bien las distancias evolutivas estimadas con las dataciones?
Para ello nos fijamos en el arbol NJ o UPGMA obtenido y en la "longitud de las 
ramas". 
Las secuencias mas antiguas deberian estar mas cerca de los nodos comunes que 
las d elos humanos, ya que tuvieron menos tiempo para acumular mutaciones. 
#¿Qué secuencias crees que se ajustan menos a su supuesta edad?
la conservacion del DNA a veces no es del todo fiable, es por eso que puede que
el material genetico haya sufrido daños que el sistema clasifica como mutaciones
estas mutaciones se pueden interpretar como ramas mas largas de lo que deberían. 
#¿Qué razones puede haber para que la edad de una secuencia no se 
#corresponda con la cantidad de cambios que ha acumulado?
Se podría deber a la misma explicacion que hemos dado en la pregunta anterior.

#Cantidad de cambios cromosomicos entre humanos y neandertales: 
distMat <- as.matrix(distancias)
# Esto son vectores lógicos que podemos usar después para seleccionar filas
# y columnas de la matriz de distancias.
modernos     <- startsWith(colnames(distMat), 'H.s.modern')
neandertales <- startsWith(colnames(distMat), 'H.s.neandertal')
denisovanos  <- startsWith(colnames(distMat), 'H.s.denisova')
chimpances   <- startsWith(colnames(distMat), 'P.troglodytes')
bonobos      <- startsWith(colnames(distMat), 'P.paniscus')

#distancia media: 
mean(distMat[modernos, neandertales])

#EJERCICIO: 
#Calcula las distancias medias entre todos los grupos.
# Distancia Modernos - Neandertales
dist_mod_nean <- mean(distMat[modernos, neandertales])

# Distancia Modernos - Denisovanos
dist_mod_deni <- mean(distMat[modernos, denisovanos])

# Distancia Modernos - Chimpancés
dist_mod_chimp <- mean(distMat[modernos, chimpances])

# Distancia Neandertales - Denisovanos
dist_nean_deni <- mean(distMat[neandertales, denisovanos])
#¿Cuántas veces mayor es la distancia entre humanos y chimpancés que entre 
#humanos y neandertales?
Dividimos la distancia humano-chimpance entre la de humano-neandertal
ratio <- dist_mod_chimp / dist_mod_nean
ratio
El resultado nos indica que un chimpancé esta 7.5 veces mas alejado de los 
humanos que un neandertal. 
#¿Cómo calcularías la distancia media entre dos secuencias de humanos modernos?
z <- distMat[modernos, modernos]
z[z == 0] <- NA
mean(z, na.rm = TRUE)
#En los árboles las distancias se representan mediante la suma de las longitudes 
#que separan las hojas. Las podemos obtener con la función:
distCofen <- cophenetic(mt_nj_rooted)

#EJERCICIO: 
#¿cuanto se parecen las distancias estimadas a las distancias cofenéticas?
# Aquí usamos las matricees como vectores:
plot(distMat, distCofen, xlab = 'Distancias F81', ylab = 'Dist. cofenéticas')
abline(a = 0, b = 1, col = 'red', lty = 2)

#limpieza o filtrado: 
crudas <- dist.dna(as.DNAbin(mtdna), model = 'raw')
crudas <- as.matrix(crudas)
identicas <- which(crudas == 0 & lower.tri(crudas), arr.ind = TRUE)
# posiciones de la semimatriz inferior donde las distancias son 0:
identicas
repetidas <- unique(c(identicas[, 'row'], identicas[, 'col']))
heatmap(crudas[repetidas, repetidas], cexRow = 0.7, cexCol = 0.7)
# La exclamación ("!") invierte el valor lógico de lo que la sigue.
# La función "grepl()" produce un vector lógico con "TRUE" donde los nombres
# de las secuencias contienen alguna de esas palabras.
mtdna37 <- mtdna[! grepl('(KX198085|KX198086|KX198082|KX198088|KX198083)',
                         names(mtdna)), ]
mtdna37

#EJERCICIO: 
#¿serias capaz de crear otro conjunto de datos, mtdna28, que solo contengan 
#secuencias del genero homo?
mtdna28 <- mtdna37[grepl('^H.', names(mtdna37)), ]
mtdna36 <- mtdna37[! startsWith(names(mtdna37), 'G.gorilla'),]
