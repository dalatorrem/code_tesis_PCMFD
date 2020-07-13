# correr si se tienen instaladas las librerias de funciones_tesis_version_2
# si no se tienen instaladas las librerias se deben instalar antes
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script.dir)
source("funciones_tesis_version_2.R")


# estos datos corresponden a los descargados de la página 
# https://archive.ics.uci.edu/ml/machine-learning-databases/00274/
# https://archive.ics.uci.edu/ml/datasets/SML2010
# data_1 corresponde a los descargado de los datos 1 analogamente para 2

# abriendo la carpeta que contiene los datos descargados de la pagina 
# se puede descargar la carpeta comprimida de la página y descomprimir ubicando la carpeta 
# NEW-DATA en la ubicación de este script
# !! antes de importar los datos es necesario 
# borrar los espacios iniciales que se tienen en el txt,
# si se importan así las columnas quedan movidas
# en ambos archivos se debe borrar "# y los dos espacios" o sea se debe dejar el 1 como primer 
# caracter


setwd(paste0(script.dir,"/NEW-DATA"))
data_1<-read.delim("NEW-DATA-1.T15.txt",header=TRUE,sep = " ")
data_2<-read.delim("NEW-DATA-2.T15.txt",header=TRUE,sep = " ")
setwd(script.dir)
str(data_1)
str(data_2)



################################################################################
################################################################################
################################################################################
## revisión de los datos
################################################################################
################################################################################
################################################################################

data_1$X1.Date<-as.factor(data_1$X1.Date)
data_2$X1.Date<-as.factor(data_2$X1.Date)
levels(data_1$X1.Date)
levels(data_2$X1.Date)
#### revisando cuántas mediciones hay por día
medidas_por_dia<-function(datos)
{
  longitudes<-c()
  for (fecha in levels(datos$X1.Date))
  {
    print(fecha)
    filtro<-datos[datos$X1.Date==fecha,]
    lon<-length(filtro$X1.Date)
    longitudes<-c(longitudes,lon)
    plot(longitudes)
    
  }
  print(max(longitudes))
  print(min(longitudes))
  print(length(longitudes))
  return(longitudes)
}  
medidas_por_dia(data_1)
medidas_por_dia(data_2)
unique(data_1$X1.Date)
unique(data_2$X1.Date)
### detectando las fechas que no tienen 96 mediciones
detectando_fechas<-function(datos)
{
  for (fecha in levels(datos$X1.Date))
  {
    
    filtro<-datos[datos$X1.Date==fecha,]
    lon<-length(filtro$X1.Date)
    if(lon!=96)
    {
      print(fecha)
    }
    
  }
}
detectando_fechas(data_1)
data_1<-filter(data_1,!(X1.Date=="11/04/2012"|X1.Date=="13/03/2012"))
data_1<-droplevels(data_1)         
detectando_fechas(data_1)



detectando_fechas(data_2)
data_2<-filter(data_2,!(X1.Date=="02/05/2012"))
data_2<-droplevels(data_2)         
detectando_fechas(data_1)
detectando_fechas(data_2)
dim(data_1)[1]/96
dim(data_2)[1]/96
# todos los registros diarios quedaron de 96 dias 
##############################################################################
##### revision de horas
##############################################################################

### revision data_1
horas<-matrix(data_1$X2.Time,nrow =96)
n<-dim(horas)[2]
m<-dim(horas)[1]
for (i in 1:m)
{
  hora<-horas[i,]
  #Sys.sleep(0.2)
  hora_revisar<-hora[1]
  print(hora_revisar)
  for (j in 1:n) 
  {
    if(hora_revisar!=hora[j])
    {
      print(c(i,j))
      print(hora[j])
      break
    }    
  }
}
### revision data_2
horas<-matrix(data_2$X2.Time,nrow =96)
n<-dim(horas)[2]
m<-dim(horas)[1]
for (i in 1:m)
{
  hora<-horas[i,]
  #Sys.sleep(0.2)
  hora_revisar<-hora[1]
  print(hora_revisar)
  for (j in 1:n) 
  {
    if(hora_revisar!=hora[j])
    {
      print(c(i,j))
      print(hora[j])
      break
    }    
  }
}
### las horas son las mismas, el for no se rompió en ninguno de los dos casos



##############################################################################
### selección de variables
##############################################################################
### data_1
str(data_1)
for (name in colnames(data_1)) 
{
  plot(data_1[,name],main=name)
  #Sys.sleep(0.2)
}
# posibes variables para el análisis
### 23 humedad exterior
### 22 temperatura exterior
### 4 temperatura habitacion
### 3 temperatura comedor sensor

### data_2
str(data_2)
for (name in colnames(data_2)) 
{
  plot(data_2[,name],main=name)
  #Sys.sleep(0.2)
}
# posibes variables para el análisis
### 23 humedad exterior
### 22 temperatura exterior
### 4 temperatura habitacion
### 3 temperatura comedor sensor
### 5 weather temperature

#### pueden ser consideradas las variables 23,22,4,3y 5 
### sin embargo en los datos 1 la variables 5 tiene un comportamiento extraño


#############################################################################

#############################################################################





################################################################################
################################################################################
################################################################################
### INICIO ANÁLISIS  # Temperaturas exterior y habitación
################################################################################
################################################################################
################################################################################

################################################################################
##################### pruebas
## los primeros 7 datos de la base de datos 1 son eliminados ya que presentan algunas
## incosistencias sobre todo el día 7, por esta razón se toman 25 dias consecutivos de 
## la primera baase y 14 consecutivos de la segunda
## hay 35 datos los primeros 21 corresponden a los datos de la base 1, los otros
## 14 a la base 2, están ordenados cronologicamente, el cambio debería estar
## en el dato 21
################################################################################




################################################################################
################### punto de cambio univariado 
################################################################################

t<-(0:95)/96


#### Habitación

temp_hab_1<-matrix(data_1$X4.Temperature_Habitacion_Sensor,nrow =96)
temp_hab_2<-matrix(data_2$X4.Temperature_Habitacion_Sensor,nrow =96)
temp_hab<-cbind(temp_hab_1,temp_hab_2)
temperatura_hab_para_pruebas<-temp_hab[,8:dim(temp_hab)[2]]

### Exterior
temperatura_1<-matrix(data_1$X22.Temperature_Exterior_Sensor,nrow =96)
temperatura_2<-matrix(data_2$X22.Temperature_Exterior_Sensor,nrow =96)
temperatura<-cbind(temperatura_1,temperatura_2)
temperatura_exterior_para_pruebas<-temperatura[,8:dim(temperatura)[2]]


# creando csv para app

#datos<-data.frame("Temp Hab"=matrix(temperatura_hab_para_pruebas,ncol = 1),
#                  "Temp Ext"=matrix(temperatura_exterior_para_pruebas,ncol = 1))
#write.table(datos,'datos_casa_domotica.csv',sep = ';',row.names = FALSE)






# plots
layout(matrix(c(1,2),nrow=1))
plot(temperatura_hab_para_pruebas[,1],type="l",col=1,
     ylab="Temperatura °C",xlab="t",
     ylim=c(min(temperatura_exterior_para_pruebas),max(temperatura_exterior_para_pruebas)))
for (j in 2:dim(temperatura_exterior_para_pruebas)[2])
{
  lines(temperatura_hab_para_pruebas[,j],col=j)
}






plot(x=t,temperatura_exterior_para_pruebas[,1],type="l",col=1,
     ylab="Temperatura °C",xlab="t",
     ylim=c(min(temperatura_exterior_para_pruebas),
            max(temperatura_exterior_para_pruebas)))
for (j in 2:dim(temperatura_exterior_para_pruebas)[2])
{
  lines(x=t,temperatura_exterior_para_pruebas[,j],col=j)
}



no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
d_no_comp_prin<-2



######## análisis habitación



temperatura_habitacion_suavizada<-Data2fd(argvals = t,temperatura_hab_para_pruebas,
                                          basisobj=Bfourier,dfscale = "gcv")

plot(temperatura_habitacion_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylim=c(0,40),ylab=c("Temperatura °C"))
comp_principales_temperatura_habitacion<-pca.fd(temperatura_habitacion_suavizada,
                                                nharm = d_no_comp_prin,
                                                harmfdPar = fdPar(Bfourier),
                                                centerfns = TRUE)
sum(comp_principales_temperatura_habitacion$varprop)
est_s_hab<-estadistico_S_N_D(comp_principales_temperatura_habitacion$scores,
                         comp_principales_temperatura_habitacion$values)
estT_hab<-estadistico_T_N_x(comp_principales_temperatura_habitacion$scores,
                            comp_principales_temperatura_habitacion$values,FALSE)


max(estT_hab)
which.max(estT_hab)
est_s_hab
sum(estT_hab)/dim(temperatura_hab_para_pruebas)[2]



######## análisis exterior

temperatura_extrerior_suavizada<-Data2fd(argvals = t,temperatura_exterior_para_pruebas,
                       basisobj=Bfourier,dfscale = "gcv")
plot(temperatura_extrerior_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylim=c(0,40),ylab=c("Temperatura °C"))
comp_principales_temperatura_externa<-pca.fd(temperatura_extrerior_suavizada,
                                          nharm = d_no_comp_prin,
                                          harmfdPar = fdPar(Bfourier),
                                          centerfns = TRUE)
sum(comp_principales_temperatura_externa$varprop)
est_s_ext<-estadistico_S_N_D(comp_principales_temperatura_externa$scores,
                  comp_principales_temperatura_externa$values)
estT_ext<-estadistico_T_N_x(comp_principales_temperatura_externa$scores,
                        comp_principales_temperatura_externa$values,TRUE)
max(estT_ext)
which.max(estT_ext)
est_s_ext
sum(estT_ext)/dim(temperatura_exterior_para_pruebas)[2]




plot(x=(0:(dim(estT_hab)[1]-1))/dim(estT_hab)[1],estT_hab,type="l",xlab="",ylab="")
plot(x=(0:(dim(estT_ext)[1]-1))/dim(estT_ext)[1],estT_ext,type="l",xlab="",ylab="")



##############################################################################################
##############################################################################################
############ metodología bivariada   
##############################################################################################
##############################################################################################
# se realiza un cambio de base para los errores, para quitar la correlación entre las variables
# el primer parametro es un arreglo de muestras de procesos bivariados
# las tres dimensiones del arreglo son:
# 1. no. de intervalos en el tiempo
# 2. dos (este es una análisis de dos procesos)
# 3. el tamaño de la muestra (la cantidad de procesos bivariados)
##############################################################################################
n<-35
no_tiempos<-96
dim(temperatura_exterior_para_pruebas) #verificando dimensiones :)
dim(temperatura_hab_para_pruebas)      #verificando dimensiones :)
# generando arreglo 
# proceso1 temperatura exterior
# proceso2 temperatura hab
procesos<-array((1:(2*n*no_tiempos)),c(no_tiempos,2,n))
for (j in 1:n)
{
  procesos[,1,j]<-temperatura_exterior_para_pruebas[,j]
  procesos[,2,j]<-temperatura_hab_para_pruebas[,j]
}



lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
##### gráficos de la transformación, valores de la matriz de covarianza



raiz_val_propios<-sqrt(lista_cambio_de_base$array_valores_pro)  

# raiz de los valores propios en el tiempo
ymin<-min(raiz_val_propios)
ymax<-max(raiz_val_propios)
plot(t,raiz_val_propios[1,],"l",col="gray",ylim = c(ymin,ymax))
lines(t,raiz_val_propios[2,],"l",col="blue")
proceso_transformado<-lista_cambio_de_base$nuevo_proceso




no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Componente 1 suavizamientos")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Componente 2 suavizamientos")


d1<-2
d2<-5

comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

sum(comp_principales_1$varprop)
sum(comp_principales_2$varprop)

T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#



T_estadistico<-T_estadistico_1+T_estadistico_2
plot(seq(1,length(T_estadistico[,1]),1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico)
which.max(T_estadistico)
(sum(T_estadistico)/35)

#########################################################################################
##### primer cambio en 21, busqueda de otros cambios:
#########################################################################################



#### análisis del dato 1 al 21.


n<-21
no_tiempos<-96

# generando arreglo 
# proceso1 temperatura exterior
# proceso2 temperatura hab
procesos<-array((1:(2*n*no_tiempos)),c(no_tiempos,2,n))
for (j in 1:21)
{
  procesos[,1,j]<-temperatura_exterior_para_pruebas[,j]
  procesos[,2,j]<-temperatura_hab_para_pruebas[,j]
}



lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
##### gráficos de la transformación, valores de la matriz de covarianza



raiz_val_propios<-sqrt(lista_cambio_de_base$array_valores_pro)  

# raiz de los valores propios en el tiempo
ymin<-min(raiz_val_propios)
ymax<-max(raiz_val_propios)
plot(t,raiz_val_propios[1,],"l",col="gray",ylim = c(ymin,ymax))
lines(t,raiz_val_propios[2,],"l",col="blue")
proceso_transformado<-lista_cambio_de_base$nuevo_proceso




no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Componente 1 suavizamientos")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Componente 2 suavizamientos")


d1<-2
d2<-5

comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

sum(comp_principales_1$varprop)
sum(comp_principales_2$varprop)

T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#



T_estadistico<-T_estadistico_1+T_estadistico_2
plot(seq(1,length(T_estadistico[,1]),1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico)
which.max(T_estadistico)
(sum(T_estadistico)/21)
### no se rechaza




#### análisis del dato 21 al 35.


n<-14
no_tiempos<-96

# generando arreglo 
# proceso1 temperatura exterior
# proceso2 temperatura hab
procesos<-array((1:(2*n*no_tiempos)),c(no_tiempos,2,n))
for (j in 22:35)
{
  procesos[,1,j-21]<-temperatura_exterior_para_pruebas[,j]
  procesos[,2,j-21]<-temperatura_hab_para_pruebas[,j]
}


lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
##### gráficos de la transformación, valores de la matriz de covarianza



raiz_val_propios<-sqrt(lista_cambio_de_base$array_valores_pro)  

# raiz de los valores propios en el tiempo
ymin<-min(raiz_val_propios)
ymax<-max(raiz_val_propios)
plot(t,raiz_val_propios[1,],"l",col="gray",ylim = c(ymin,ymax))
lines(t,raiz_val_propios[2,],"l",col="blue")
proceso_transformado<-lista_cambio_de_base$nuevo_proceso




no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Componente 1 suavizamientos")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Componente 2 suavizamientos")


d1<-2
d2<-5

comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

sum(comp_principales_1$varprop)
sum(comp_principales_2$varprop)

T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#



T_estadistico<-T_estadistico_1+T_estadistico_2
plot(seq(1,length(T_estadistico[,1]),1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico)
which.max(T_estadistico)
(sum(T_estadistico)/14)
### no se rechaza



