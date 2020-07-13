### simulación de proceso p-variantes
direccion_script<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(direccion_script)
source("funciones_tesis_version_2.R")
# definir los puntos muestreados en el intervalo [0,1]
t<-(0:100)/100
p<-3
n_len_t<-length(t)
N_tamano<-50
####################################################################################################
# defininiendo p^2 funciones para las componentes de la matriz de transformación en este caso p=3
####################################################################################################
p1<-t*2
p2<-t+1            # cambie esta función de la anterior simulación
p3<-t^2            # cambie esta función de la anterior simulación
p4<-exp(t)         
p5<-sin(3*t)+2     # cambie esta función de la anterior simulación
p6<-cos(4+t)+2     # cambie esta función de la anterior simulación
p7<-t^3            # cambie esta función de la anterior simulación
p8<-0*t+1.5        # cambie esta función de la anterior simulación
p9<-t^2+1          # cambie esta función de la anterior simulación
####################################################################################################
## Se colocan las  p^2 componentes dentro de una matriz
####################################################################################################
mat_componentes<-matrix(c(p1,p2,p3,p4,p5,p6,p7,p8,p9),nrow=n_len_t)

####################################################################################################
## se verifica que los valores propios de A*A' no son cero y que las funciones
## de  los valores propios de estos no se cruzan
####################################################################################################
arreglo_de_transformaciones<-verif_arreg_de_transformaciones(mat_componentes,t)

####################################################################################################
#  se incorporan p-veces el argumento t
####################################################################################################
argvalsList <- list( t, t,t)

####################################################################################################
# se fija la semilla y se simula con la función del artículo de clara happ 
####################################################################################################
set.seed(1234)
objeto_de_Y_sin_ruido <- simMultiFunData(N = N_tamano, argvals = argvalsList,
                                         eFunType = "Fourier", eValType = "linear", M = 10,
                                         type = "split")$simData
# Extrayendo el valor a mayor distancia del 0 de los simulados de Y 
max_dist_Y_sin_ruido<-max(sqrt(rbind(objeto_de_Y_sin_ruido@.Data[[1]]@X,
                                     objeto_de_Y_sin_ruido@.Data[[2]]@X,
                                     objeto_de_Y_sin_ruido@.Data[[3]]@X)^2))
max_dist_Y_sin_ruido  #  4.823947
## fijando como desviacion el max extraido
desv<-max_dist_Y_sin_ruido
# fijando semilla para los ruidos
set.seed(1234)
arreglo_de_Y_con_ruido<-addError(objeto_de_Y_sin_ruido,sd=desv)
plot(arreglo_de_Y_con_ruido,xlab="")

####################################################################################################
### se deja que el maximo error sea 1 (solamente por convención)
####################################################################################################
max_dist_Y_con_ruido<-max(sqrt(rbind(arreglo_de_Y_con_ruido@.Data[[1]]@X,
                                     arreglo_de_Y_con_ruido@.Data[[2]]@X,
                                     arreglo_de_Y_con_ruido@.Data[[3]]@X)^2))
max_dist_Y_con_ruido 
arreglo_de_Y_con_ruido@.Data[[1]]@X<-arreglo_de_Y_con_ruido@.Data[[1]]@X/max_dist_Y_con_ruido
arreglo_de_Y_con_ruido@.Data[[2]]@X<-arreglo_de_Y_con_ruido@.Data[[2]]@X/max_dist_Y_con_ruido
arreglo_de_Y_con_ruido@.Data[[3]]@X<-arreglo_de_Y_con_ruido@.Data[[3]]@X/max_dist_Y_con_ruido
max_dist_Y_con_ruido<-max(sqrt(rbind(arreglo_de_Y_con_ruido@.Data[[1]]@X,
                                     arreglo_de_Y_con_ruido@.Data[[2]]@X,
                                     arreglo_de_Y_con_ruido@.Data[[3]]@X)^2))
max_dist_Y_con_ruido #1
plot(arreglo_de_Y_con_ruido,xlab="")
####################################################################################################
## boxplot para los tres ruidos
####################################################################################################
matrix_para_boxplot<-matrix(arreglo_de_Y_con_ruido@.Data[[1]]@X,ncol = 1)
matrix_para_boxplot<-cbind(matrix_para_boxplot,matrix(arreglo_de_Y_con_ruido@.Data[[2]]@X,ncol = 1))
matrix_para_boxplot<-cbind(matrix_para_boxplot,matrix(arreglo_de_Y_con_ruido@.Data[[3]]@X,ncol = 1))
sd(arreglo_de_Y_con_ruido@.Data[[1]]@X)
sd(arreglo_de_Y_con_ruido@.Data[[2]]@X) 
sd(arreglo_de_Y_con_ruido@.Data[[3]]@X) 
boxplot(matrix_para_boxplot)
####################################################################################################
## Función de media antes del cambio
####################################################################################################
fun_media_1_parte_1<-2*(t-0.5)^2
fun_media_2_parte_1<-3*t^3-0.9*t
fun_media_3_parte_1<-log(t+1)


####################################################################################################
# Función de medias después primer  cambio
####################################################################################################
v1<-0.05 #
v2<-0.05   # cambie esta función de la simulación anterior
v3<-0.05   # cambie esta función de la simulación anterior
# 
fun_media_1_parte_2<-2*(t-0.5)^2+v1*t^2
fun_media_2_parte_2<-3*t^3-0.9*t+v2*exp(t)
fun_media_3_parte_2<-log(t+1)+v3


####################################################################################################
# Función de medias después segundo  cambio
####################################################################################################
v1<-0.05 #
v2<-0.05   # cambie esta función de la simulación anterior
v3<-0.05  # cambie esta función de la simulación anterior
# 
fun_media_1_parte_3<-2*(t-0.5)^2+v1*cos(8*pi*t)
fun_media_2_parte_3<-3*t^3-0.9*t+v2*t
fun_media_3_parte_3<-(1+v3)*log(t+1)

Conf3x2 = matrix(c(1:3), nrow=1, byrow=TRUE)
layout(Conf3x2)


y_min=min(fun_media_1_parte_1,fun_media_1_parte_2,fun_media_1_parte_3)
y_max=max(fun_media_1_parte_1,fun_media_1_parte_2,fun_media_1_parte_3)
plot(t,fun_media_1_parte_1,"l",col=rgb(0,0,0),main="",ylab="",ylim=c(y_min,y_max))
lines(t,fun_media_1_parte_2,"l",col=rgb(0.5,0.5,0.5))
lines(t,fun_media_1_parte_3,"l",col=rgb(0,0,1))

y_min=min(fun_media_2_parte_1,fun_media_2_parte_2,fun_media_2_parte_3)
y_max=max(fun_media_2_parte_1,fun_media_2_parte_2,fun_media_2_parte_3)
plot(t,fun_media_2_parte_1,"l",col=rgb(0,0,0),main="",ylab="",ylim=c(y_min,y_max))
lines(t,fun_media_2_parte_2,"l",col=rgb(0.5,0.5,0.5))
lines(t,fun_media_2_parte_3,"l",col=rgb(0,0,1))

y_min=min(fun_media_3_parte_1,fun_media_3_parte_2,fun_media_3_parte_3)
y_max=max(fun_media_3_parte_1,fun_media_3_parte_2,fun_media_3_parte_3)
plot(t,fun_media_3_parte_1,"l",col=rgb(0,0,0),main="",ylab="",ylim=c(y_min,y_max))
lines(t,fun_media_3_parte_2,"l",col=rgb(0.5,0.5,0.5))
lines(t,fun_media_3_parte_3,"l",col=rgb(0,0,1))

# punto de cambio debe ser menor a N_tamano
n1<-15
n2<-24
n3<-N_tamano-n1-n2
####################################################################################################
# arreglos con las funciones medias
####################################################################################################
arreglo_de_fun_media_1_parte_1<-kronecker(t(rep(1,n1)),fun_media_1_parte_1)
arreglo_de_fun_media_2_parte_1<-kronecker(t(rep(1,n1)),fun_media_2_parte_1)
arreglo_de_fun_media_3_parte_1<-kronecker(t(rep(1,n1)),fun_media_3_parte_1)
##
arreglo_de_fun_media_1_parte_2<-kronecker(t(rep(1,n2)),fun_media_1_parte_2)
arreglo_de_fun_media_2_parte_2<-kronecker(t(rep(1,n2)),fun_media_2_parte_2)
arreglo_de_fun_media_3_parte_2<-kronecker(t(rep(1,n2)),fun_media_3_parte_2)
##
arreglo_de_fun_media_1_parte_3<-kronecker(t(rep(1,n3)),fun_media_1_parte_3)
arreglo_de_fun_media_2_parte_3<-kronecker(t(rep(1,n3)),fun_media_2_parte_3)
arreglo_de_fun_media_3_parte_3<-kronecker(t(rep(1,n3)),fun_media_3_parte_3)
##
arreglo_de_fun_media_1<-cbind(arreglo_de_fun_media_1_parte_1,arreglo_de_fun_media_1_parte_2,
                              arreglo_de_fun_media_1_parte_3)
arreglo_de_fun_media_2<-cbind(arreglo_de_fun_media_2_parte_1,arreglo_de_fun_media_2_parte_2,
                              arreglo_de_fun_media_2_parte_3)
arreglo_de_fun_media_3<-cbind(arreglo_de_fun_media_3_parte_1,arreglo_de_fun_media_3_parte_2,
                              arreglo_de_fun_media_3_parte_3)
####################################################################################################
## generando arreglo de los Y independientes
#################################################################################################
arreglo_Y_gorro<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
arreglo_Y_gorro[,1,]<-t(arreglo_de_Y_con_ruido@.Data[[1]]@X)
arreglo_Y_gorro[,2,]<-t(arreglo_de_Y_con_ruido@.Data[[2]]@X)
arreglo_Y_gorro[,3,]<-t(arreglo_de_Y_con_ruido@.Data[[3]]@X)

#################################################################################################
#### transformando los errores ( se toma arreglo de transformaciones que sale de la función )
#################################################################################################
max_transformacion<-sqrt(max(arreglo_de_transformaciones^2))
arreglo_Y<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
for (j in 1:n_len_t)
{
  arreglo_Y[j,,]<-arreglo_de_transformaciones[,,j]%*%arreglo_Y_gorro[j,,]
}
#################################################################################################
## dividiendo el arreglo entre la máxima entrada para que todas las entradas 
## queden con errores entre -1 y 1
#################################################################################################
arreglo_Y<-arreglo_Y/max_transformacion
datos<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
datos[,1,]<-arreglo_de_fun_media_1 + arreglo_Y[,1,]
datos[,2,]<-arreglo_de_fun_media_2 + arreglo_Y[,2,]
datos[,3,]<-arreglo_de_fun_media_3 + arreglo_Y[,3,]

#################################################################################################
## gráficos errores
#################################################################################################
# gráfico proceso 1 con dependencia
plot(t,arreglo_Y[,1,1],"l",col=1,main="",ylab="",ylim=c(-1,1))# errores proceso 1 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,1,i],"l",col=gray(1-0.8*(i/N_tamano)))
  if(i==n1)
  {
    #Sys.sleep(1)
  }
}
# gráfico proceso 2 con dependencia
plot(t,arreglo_Y[,2,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#errores proceso 2 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,2,i],"l",col=gray(1-0.8*(i/N_tamano)))
}
# gráfico proceso 3 con dependencia
plot(t,arreglo_Y[,3,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#errores proceso 3 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,3,i],"l",col=gray(1-0.8*(i/N_tamano)))
}


#################################################################################################
## gráficos procesos 
#################################################################################################
# gráfico proceso 1 con dependencia
y_max=max(datos)
y_min=min(datos)
plot(t,datos[,1,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X1",ylim=c(y_min,y_max))
for (i in 2:N_tamano)
{
  escala_color=0.8*(N_tamano-i)/N_tamano+0.2
  lines(t,datos[,1,i],"l",col=rgb(escala_color,escala_color,1))
}
plot(t,datos[,2,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X2",ylim=c(y_min,y_max))
for (i in 2:N_tamano)
{
  escala_color=0.8*(N_tamano-i)/N_tamano+0.2
  lines(t,datos[,2,i],"l",col=rgb(escala_color,escala_color,1))
}
plot(t,datos[,3,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X3",ylim=c(y_min,y_max))
for (i in 2:N_tamano)
{
  escala_color=0.8*(N_tamano-i)/N_tamano+0.2
  lines(t,datos[,3,i],"l",col=rgb(escala_color,escala_color,1))
}
#################################################################################################
#################################################################################################
#### TERMINA LA SIMULACIÓN DE LOS DATOS
#################################################################################################
#################################################################################################
#### EMPIEZA EL ANÁLISIS
#################################################################################################
#################################################################################################


#################################################################################################
### encontrando componentes principales
#################################################################################################
lista_com_principales<-componentes_principales_multivariadas_continuas(datos)
lista_com_principales$array_mat_covarianzas
# valores propios
valores_propios<-lista_com_principales$array_valores_pro
ymax<-max(log(valores_propios))
ymin<-min(log(valores_propios))
layout(1)
verif_arreg_de_transformaciones(mat_componentes,t)
plot(t,log(valores_propios[1,]),"p",col=1,
     main="",ylim=c(ymin,ymax),ylab="")#Valores propios en el tiempo
lines(t,log(valores_propios[2,]),"p",col=2)
lines(t,log(valores_propios[3,]),"p",col=3)
ajuste_val_1<-lm(log(valores_propios[1,])~bs(t,5))
ajuste_val_2<-lm(log(valores_propios[2,])~bs(t,5))
ajuste_val_3<-lm(log(valores_propios[3,])~bs(t,5))

lines(t,ajuste_val_1$fitted.values,col=1)
lines(t,ajuste_val_2$fitted.values,col=2)
lines(t,ajuste_val_3$fitted.values,col=3)
#################################################################################################
# haciendo el cambio de base que hace las matrices de covarianza sean identidades
#################################################################################################
lista_transformados<-cambio_de_base_comp_principales(datos)
medias<-lista_transformados$arreglo_vec_medias

#################################################################################################
## gráficos media teórica y media estimada
#################################################################################################

plot(t,fun_media_1_parte_1,"l",col="green",main="Medias proceso 1")
lines(t,fun_media_1_parte_2,"l",col="blue")
lines(t,medias[,1],"l",col="red")
lines(t,fun_media_1_parte_3,"l",col="yellow")
plot(t,fun_media_2_parte_1,"l",col="green",main="Medias proceso 2")
lines(t,fun_media_2_parte_2,"l",col="blue")
lines(t,medias[,2],"l",col="red")
lines(t,fun_media_2_parte_3,"l",col="yellow")
plot(t,fun_media_3_parte_1,"l",col="green",main="Medias proceso 2")
lines(t,fun_media_3_parte_2,"l",col="blue")
lines(t,medias[,3],"l",col="red")
lines(t,fun_media_3_parte_3,"l",col="yellow")


#################################################################################################
# proceso transformado
#################################################################################################
proceso_transformado<-lista_transformados$nuevo_proceso
# las matrices de covarianza en cada punto deberían ser la identidad y los valores propios todos 
# deberían ser 1

#########################################################################################################
## suavizando cada una de las componentes
#########################################################################################################
no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Datos suavizados")
#########################################################################################################
### análisis de punto de cambio
#########################################################################################################
# d no. de componentes principales
d1<-10
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#

d2<-10
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
d3<-8

comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)
plot(seq(1,N_tamano,1)/N_tamano,T_estadistico,main="",ylab="","l",xlab="")
max(T_estadistico)
which.max(T_estadistico)




sum(T_estadistico)/N_tamano
which.max(T_estadistico)  





#################################################################################################
#################################################################################################
### Busqueda de otros cambios
#################################################################################################
#################################################################################################

#################################################################################################
# parte_1
#################################################################################################
dim(datos)
datos_1<-datos[,,1:39]


lista_transformados<-cambio_de_base_comp_principales(datos_1)
proceso_transformado<-lista_transformados$nuevo_proceso

no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Datos suavizados")
# d no. de componentes principales
d1<-10
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#

d2<-10

comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
d3<-8
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#


T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)

max(T_estadistico)
which.max(T_estadistico)


#

sum(T_estadistico)/length(T_estadistico)
which.max(T_estadistico) # se rechaza hipotesis nula 
plot(seq(1,length(T_estadistico),1),T_estadistico,main="",ylab="","l",xlab="")


#################################################################################################
# parte_1_1
#################################################################################################
dim(datos)
datos_1_1<-datos[,,1:15]


lista_transformados<-cambio_de_base_comp_principales(datos_1_1)
proceso_transformado<-lista_transformados$nuevo_proceso

no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Datos suavizados")
# d no. de componentes principales
d1<-8
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#
d2<-9
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
d3<-7


comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)

max(T_estadistico)
which.max(T_estadistico)


#

sum(T_estadistico)/16
which.max(T_estadistico) # NO se rechaza hipotesis nula 
plot(seq(1,length(T_estadistico),1),T_estadistico,main="",ylab="","l",xlab="")
#################################################################################################
# parte_1_2
#################################################################################################
dim(datos)
datos_1_2<-datos[,,16:39]


lista_transformados<-cambio_de_base_comp_principales(datos_1_2)
proceso_transformado<-lista_transformados$nuevo_proceso

no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Datos suavizados")
# d no. de componentes principales
d1<-9
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#
d2<-10
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
d3<-9


comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)

max(T_estadistico)
which.max(T_estadistico)


#

sum(T_estadistico)/23
which.max(T_estadistico) # NO se rechaza hipotesis nula 

plot(seq(1,length(T_estadistico),1),T_estadistico,main="",ylab="","l",xlab="")


#################################################################################################
# parte_2
#################################################################################################
dim(datos)
datos_2<-datos[,,40:50]


lista_transformados<-cambio_de_base_comp_principales(datos_2)
proceso_transformado<-lista_transformados$nuevo_proceso

no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Datos suavizados")
# d no. de componentes principales
d1<-7
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#


d2<-6
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#

d3<-7


comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)

max(T_estadistico)
which.max(T_estadistico)


#
sum(T_estadistico)/11
which.max(T_estadistico) # NO se rechaza hipotesis nula 
plot(seq(1,length(T_estadistico),1),T_estadistico,main="",ylab="","l",xlab="")
