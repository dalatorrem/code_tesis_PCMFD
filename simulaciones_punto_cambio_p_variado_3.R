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
fun_media_1_antes<-t +exp(2*t)
fun_media_2_antes<-2*cos(3*(2*pi)*t)+sin(9*(2*pi)*t)
fun_media_3_antes<-2*sin((2*pi)*t)+sin(5*(2*pi)*t)+-1


layout(matrix(c(1,2,3),nrow=1))
layout.show(3)

plot(t,fun_media_1_antes,type="l",col=1,ylab="")
plot(t,fun_media_2_antes,type="l",col=2,ylab="")
plot(t,fun_media_3_antes,type="l",col=3,ylab="")
####################################################################################################
# Función de medias después del cambio
####################################################################################################
v1<-0.05 #
v2<-0.05 # cambie esta función de la simulación anterior
v3<-0.05# cambie esta función de la simulación anterior
# 

### la función PRUEBA se definió para correr rapidamente todo el análisis, 
### sin emabrgo es necesario observar con cuidado todo lo que se hace al interior de ésta.
### la función retorna los datos simulados
datos<-PRUEBA(10,9,9)


PRUEBA<-function(d1,d2,d3)
{
# definiendo las funciones medias después del punto de cambio
fun_media_1_despues<-t +exp(2*t)+v1*t^2
plot(fun_media_1_despues,type="l")
fun_media_2_despues<-2*cos(3*(2*pi)*t)+(1+v2)*(sin(9*(2*pi)*t))
plot(fun_media_2_despues,type="l")
fun_media_3_despues<-2*sin((2*pi)*t)+sin(5*(2*pi)*t+v3*(2*pi*t))-1
plot(fun_media_3_despues,type="l")
# punto de cambio debe ser menor a N_tamano
n1<-15
n2<-N_tamano-n1
####################################################################################################
# arreglos con las funciones medias
####################################################################################################
arreglo_de_fun_media_1_antes<-kronecker(t(rep(1,n1)),fun_media_1_antes)
arreglo_de_fun_media_2_antes<-kronecker(t(rep(1,n1)),fun_media_2_antes)
arreglo_de_fun_media_3_antes<-kronecker(t(rep(1,n1)),fun_media_3_antes)
##
arreglo_de_fun_media_1_despues<-kronecker(t(rep(1,n2)),fun_media_1_despues)
arreglo_de_fun_media_2_despues<-kronecker(t(rep(1,n2)),fun_media_2_despues)
arreglo_de_fun_media_3_despues<-kronecker(t(rep(1,n2)),fun_media_3_despues)
##
arreglo_de_fun_media_1<-cbind(arreglo_de_fun_media_1_antes,arreglo_de_fun_media_1_despues)
arreglo_de_fun_media_2<-cbind(arreglo_de_fun_media_2_antes,arreglo_de_fun_media_2_despues)
arreglo_de_fun_media_3<-cbind(arreglo_de_fun_media_3_antes,arreglo_de_fun_media_3_despues)
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
# gráfico proceso 2 con dependencia
plot(t,arreglo_Y[,1,1],"l",col=1,main="",ylab="",ylim=c(-1,1))# gráfico proceso 1 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,1,i],"l",col=i)
  if(i==n1)
  {
    #Sys.sleep(1)
  }
}
# gráfico proceso 2 con dependencia
plot(t,arreglo_Y[,2,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#camintas proceso 2 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,2,i],"l",col=i)
}
# gráfico proceso 3 con dependencia
plot(t,arreglo_Y[,3,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#camintas proceso 3 con dependencia
for (i in 1:N_tamano)
{
  lines(t,arreglo_Y[,3,i],"l",col=i)
}


#################################################################################################
## gráficos procesos con las medias
#################################################################################################
# gráfico proceso 2 con dependencia
plot(t,datos[,1,1],"l",col=1,main="",ylab="")# gráfico proceso 1 con dependencia
for (i in 1:N_tamano)
{
  lines(t,datos[,1,i],"l",col=i)
  if(i==n1)
  {
    #Sys.sleep(1)
  }
}
# gráfico proceso 2 con dependencia
plot(t,datos[,2,1],"l",col=1,main="",ylab="")#camintas proceso 2 con dependencia
for (i in 1:N_tamano)
{
  lines(t,datos[,2,i],"l",col=i)
}
# gráfico proceso 3 con dependencia
plot(t,datos[,3,1],"l",col=1,main="",ylab="")#camintas proceso 3 con dependencia
for (i in 1:N_tamano)
{
  lines(t,datos[,3,i],"l",col=i)
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
plot(t,log(valores_propios[1,]),"l",col=1,
     main="",ylim=c(ymin,ymax),ylab="")#Valores propios en el tiempo
lines(t,log(valores_propios[2,]),"l",col=2)
lines(t,log(valores_propios[3,]),"l",col=3)
#################################################################################################
# haciendo el cambio de base que hace las matrices de covarianza sean identidades
#################################################################################################
lista_transformados<-cambio_de_base_comp_principales(datos)
medias<-lista_transformados$arreglo_vec_medias

#################################################################################################
## gráficos media teórica y media estimada
#################################################################################################

plot(t,fun_media_1_antes,"l",col="green",main="Medias proceso 1")
lines(t,fun_media_1_despues,"l",col="blue")
lines(t,medias[,1],"l",col="red")
plot(t,fun_media_2_antes,"l",col="green",main="Medias proceso 2")
lines(t,fun_media_2_despues,"l",col="blue")
lines(t,medias[,2],"l",col="red")
plot(t,fun_media_3_antes,"l",col="green",main="Medias proceso 3")
lines(t,fun_media_3_despues,"l",col="blue")
lines(t,medias[,3],"l",col="red")
## los cambios son apenas notados en las gráficas

#################################################################################################
# proceso transformado
#################################################################################################
proceso_transformado<-lista_transformados$nuevo_proceso
# las matrices de covarianza en cada punto deberían ser la identidad y los valores propios todos 
# deberían ser 1
prueba_comp<-componentes_principales_multivariadas_continuas(proceso_transformado)
prueba_comp$array_mat_covarianzas
prueba_comp$array_valores_pro
# gráficas  
layout(matrix(c(1,2,3),nrow=1))
plot(t,proceso_transformado[,1,1],"l",ylab="")
for (i in 1:N_tamano)
{
  lines(t,proceso_transformado[,1,i],col=i) 
}  
plot(t,proceso_transformado[,2,1],"l",ylab="")
for (i in 1:N_tamano)
{
  lines(t,proceso_transformado[,2,i],col=i) 
}  
plot(t,proceso_transformado[,3,1],"l",ylab="")
for (i in 1:N_tamano)
{
  lines(t,proceso_transformado[,3,i],col=i) 
}  
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
#d1<-8
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
print("Variabilidad explicada componentes de var 1,2,3; respectivamente")
print(sum(comp_principales_1$varprop))
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)


#d2<-4
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
print(sum(comp_principales_2$varprop))
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)

#d3<-11
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
print(sum(comp_principales_3$varprop))
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
layout(1)
plot(seq(1,N_tamano,1)/N_tamano,T_estadistico,main="",ylab="","l",xlab="")
max(T_estadistico)
which.max(T_estadistico)


## solo de comparación
sum(T_estadistico_1)/N_tamano
which.max(T_estadistico_1)
sum(T_estadistico_2)/N_tamano
which.max(T_estadistico_2)
sum(T_estadistico_3)/N_tamano
which.max(T_estadistico_3)


#### resultados
print("estadistico de prueba")
print(sum(T_estadistico)/N_tamano)  # estadistico de prueba para comparar con la distribucion K_d
print("punto de cambio")
print(which.max(T_estadistico))   # estimación punto de cambio
return(datos)
}
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#### comparación del uso de la propuesta de Berkes  sin hacer transformaciones,
#### es decir sin tener en cuenta la correlación
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

#########################################################################################################
## suavizando cada una de las componentes
#########################################################################################################
no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=datos[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Datos suavizados")
#
smooth_data_2<-Data2fd(argvals = t,y=datos[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Datos suavizados")
#
smooth_data_3<-Data2fd(argvals = t,y=datos[,3,],basisobj=Bfourier,dfscale = "gcv")
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
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,TRUE)


d2<-10

comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,TRUE)

d3<-9
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,TRUE)

#

#


#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
plot(seq(1,N_tamano,1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico)
which.max(T_estadistico)


#  solo de referencia
sum(T_estadistico_1)/N_tamano
which.max(T_estadistico_1)
sum(T_estadistico_2)/N_tamano
which.max(T_estadistico_2)
sum(T_estadistico_3)/N_tamano
  which.max(T_estadistico_3)

### análisis de la suma de todas los estadísticos T
sum(T_estadistico)/N_tamano
which.max(T_estadistico)

