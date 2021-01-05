#Con este script lo que hacemos es una lista de todas las combinaciones posibles entre dos intervalos
#de fechas, cada "intervalo" dias.

#Por ejemplo, si quisiesemos todas las combianciones de fechas posibles de incio y fin de un 
#periodo que empiece en febrero (cualquier dia) y acabe en abril (cualquier día), fecha_inicial y 
#fecha_final serán el primer y último día de febrero respectivamente, y fecha_incial2 y fecha_final2
#serán el primer y último día de abril respectivamente. Esto resultará en una lista de todas 
#las combinaciones de fechas posibles.

fecha_inicial <- as.Date("2020-02-01")#primer día del inicio de las fechas
fecha_final <- as.Date("2020-02-28")#último día del inicio de las fechas

fecha_inicial2 <- as.Date("2020-02-1")#primer dia del final de las fechas
fecha_final2 <- as.Date("2020-02-28")#ultimo día del final de las fechas
intervalo <- 2

int1 <- seq.Date(from = fecha_inicial,to = fecha_final, by = intervalo)
int2 <- seq.Date(from = fecha_inicial2,to = fecha_final2, by = intervalo)


a <- int1
b <- int2
final <- data.frame(matrix(ncol = 2, nrow = 0))
for ( i in 1:length(a)){
  intervals <- NULL
  intervals <- expand.grid(a =a[i] , b = b[which(b>=a[i])[1]:length(b)])
  final <- rbind(final, intervals)
  
}

#Estas columnas deberemos organizarlas en un data.frame con los nombres de columnas correspondientes
#al script base.r. El archivo "interventions_template.csv" puede servir de guia para conocer la forma
#de dicho archivo.
