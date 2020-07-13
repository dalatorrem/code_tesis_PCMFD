### simulación de proceso p-variantes
### ruidos dobl_exponencial ver linea 40
### medias que no hacen parte de la base de fourier 
direccion_script<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(direccion_script)
source("funciones_tesis_version_2.R")
# definir los puntos muestreados en el intervalo [0,1]
t<-(0:100)/100
p<-4
n_len_t<-length(t)
N_tamano<-50
# definir p^2 funciones para las componentes de la matriz de transformación en este caso p=3
p1<-t*2
p2<-t+1            
p3<-t^2            
p4<-exp(t)         
p5<-sin(3*t)+2     
p6<-cos(4+t)+2    
p7<-t^3            
p8<-0*t+1.5        
p9<-t^2+1
p10<-exp(t)/3
p11<-cos(2+t)+1
p12<-cos(4*t)+1.3
p13<-t^2+1
p14<-t^3+1
p15<-t+1
p16<-t^5+1
## colocar las p^2 componentes en el c() dentro de la matriz
mat_componentes<-matrix(c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16),nrow=n_len_t)
## verificar que los valores propios de A*A' no son cero y que las funciones de estos no se cruzan
arreglo_de_transformaciones<-verif_arreg_de_transformaciones(mat_componentes,t)
#  se incorporan p-veces el argumento t
argvalsList <- list( t, t,t,t)
# fijar semilla y simular con la función del artículo de clara happ 
set.seed(1234)
objeto_de_Y_sin_ruido <- simMultiFunData(N = N_tamano, argvals = argvalsList,
                                         eFunType = "Fourier", eValType = "linear", M = 10,
                                         type = "split")$simData
plot(objeto_de_Y_sin_ruido)
# Extrayendo el valor a mayor distancia del 0 de los simulados de Y 
max_dist_Y_sin_ruido<-max(sqrt(rbind(objeto_de_Y_sin_ruido@.Data[[1]]@X,
                                     objeto_de_Y_sin_ruido@.Data[[2]]@X,
                                     objeto_de_Y_sin_ruido@.Data[[3]]@X,
                                     objeto_de_Y_sin_ruido@.Data[[4]]@X)^2))
max_dist_Y_sin_ruido  #  4.308653
## fijando como desviacion el max extraido
desv<-max_dist_Y_sin_ruido
# fijando semilla para los ruidos con double exponential
set.seed(123)
ruidos_double_exponential<-rlaplace(n_len_t*p*N_tamano,m=0,s=desv/sqrt(2))
sd(ruidos_double_exponential) # comparar con desv!!!! ojo con el parámtero de la doble exponencial
arreglo_ruidos<-array(ruidos_double_exponential,c(N_tamano,n_len_t,p))
### arreglo de ruidos
arreglo_de_Y_con_ruido<-array(rep(0,n_len_t*p*N_tamano),c(N_tamano,n_len_t,p))
arreglo_de_Y_con_ruido[,,1]<-objeto_de_Y_sin_ruido@.Data[[1]]@X+arreglo_ruidos[,,1]
arreglo_de_Y_con_ruido[,,2]<-objeto_de_Y_sin_ruido@.Data[[2]]@X+arreglo_ruidos[,,2]
arreglo_de_Y_con_ruido[,,3]<-objeto_de_Y_sin_ruido@.Data[[3]]@X+arreglo_ruidos[,,3]
arreglo_de_Y_con_ruido[,,4]<-objeto_de_Y_sin_ruido@.Data[[4]]@X+arreglo_ruidos[,,4]
### se deja que el maximo error sea 1 solamente por convención
max_dist_Y_con_ruido<-max(sqrt((arreglo_de_Y_con_ruido)^2))
max_dist_Y_con_ruido #28.53863
arreglo_de_Y_con_ruido<-arreglo_de_Y_con_ruido/max_dist_Y_con_ruido
max_dist_Y_con_ruido<-max(sqrt((arreglo_de_Y_con_ruido)^2))
max_dist_Y_con_ruido #1

sd(arreglo_de_Y_con_ruido[,,1]) #0.1576482
sd(arreglo_de_Y_con_ruido[,,2]) #0.1543155
sd(arreglo_de_Y_con_ruido[,,3]) #0.1536216
sd(arreglo_de_Y_con_ruido[,,4]) #0.1563773


matrix_para_boxplot<-matrix(c(matrix(arreglo_de_Y_con_ruido[,,1],ncol=1),
                              matrix(arreglo_de_Y_con_ruido[,,2],ncol=1),
                              matrix(arreglo_de_Y_con_ruido[,,3],ncol=1),
                              matrix(arreglo_de_Y_con_ruido[,,3],ncol=1)),ncol=4)

boxplot(matrix_para_boxplot)

##########################################################################################################
## media
#antes
fun_media_1_antes<-6*(t-0.5)^2
fun_media_2_antes<-24*(t-0.5)^4
fun_media_3_antes<-1.5*t
fun_media_4_antes<-exp(t)-1
# cambios sin cambios 
v1<-0.025
v2<-0.025
v3<-0.025
v4<-0.025
# 
fun_media_1_despues<-6*(t-0.5)^2+v1*t
fun_media_2_despues<-(1+v2)*24*(t-0.5)^4
fun_media_3_despues<-1.5*t+v3*sin(10*pi*t)
fun_media_4_despues<-(1+v4)*(exp(t)-1)
##############################################################



  layout(matrix(c(1,2,3,4),nrow = 1))
  y_min=min(fun_media_1_antes,fun_media_2_antes,fun_media_3_antes,fun_media_4_antes,
              fun_media_1_despues,fun_media_2_despues,fun_media_3_despues,fun_media_4_despues)
  y_max=max(fun_media_1_antes,fun_media_2_antes,fun_media_3_antes,fun_media_4_antes,
            fun_media_1_despues,fun_media_2_despues,fun_media_3_despues,fun_media_4_despues)
  plot(t,fun_media_1_antes,"l",col=rgb(0.5,0.5,1),main="",ylab="",ylim=c(y_min,y_max))
  lines(t,fun_media_1_despues,"l",col=rgb(0,0,1))
  plot(t,fun_media_2_antes,"l",col=rgb(0.5,0.5,1),main="",ylab="",ylim=c(y_min,y_max))
  lines(t,fun_media_2_despues,"l",col=rgb(0,0,1))
  plot(t,fun_media_3_antes,"l",col=rgb(0.5,0.5,1),main="",ylab="",ylim=c(y_min,y_max))
  lines(t,fun_media_3_despues,"l",col=rgb(0,0,1))
  plot(t,fun_media_4_antes,"l",col=rgb(0.5,0.5,1),main="",ylab="",ylim=c(y_min,y_max))
  lines(t,fun_media_4_despues,"l",col=rgb(0,0,1))
  # punto de cambio debe ser menor a N_tamano
  n1<-40
  n2<-N_tamano-n1
  # arreglos con las funciones medias
  arreglo_de_fun_media_1_antes<-kronecker(t(rep(1,n1)),fun_media_1_antes)
  arreglo_de_fun_media_2_antes<-kronecker(t(rep(1,n1)),fun_media_2_antes)
  arreglo_de_fun_media_3_antes<-kronecker(t(rep(1,n1)),fun_media_3_antes)
  arreglo_de_fun_media_4_antes<-kronecker(t(rep(1,n1)),fun_media_4_antes)
  ##
  arreglo_de_fun_media_1_despues<-kronecker(t(rep(1,n2)),fun_media_1_despues)
  arreglo_de_fun_media_2_despues<-kronecker(t(rep(1,n2)),fun_media_2_despues)
  arreglo_de_fun_media_3_despues<-kronecker(t(rep(1,n2)),fun_media_3_despues)
  arreglo_de_fun_media_4_despues<-kronecker(t(rep(1,n2)),fun_media_4_despues)
  ##
  arreglo_de_fun_media_1<-cbind(arreglo_de_fun_media_1_antes,arreglo_de_fun_media_1_despues)
  arreglo_de_fun_media_2<-cbind(arreglo_de_fun_media_2_antes,arreglo_de_fun_media_2_despues)
  arreglo_de_fun_media_3<-cbind(arreglo_de_fun_media_3_antes,arreglo_de_fun_media_3_despues)
  arreglo_de_fun_media_4<-cbind(arreglo_de_fun_media_4_antes,arreglo_de_fun_media_4_despues)
  ##########################################################################################################
  ## caminatas  de procesos p-variadas  X, con Y_k(t) dependientes
  ## punto de cambio en n1
  # Y_gorro indica procesos independientes en cada t
  arreglo_Y_gorro<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
  arreglo_Y_gorro[,1,]<-t(arreglo_de_Y_con_ruido[,,1])
  arreglo_Y_gorro[,2,]<-t(arreglo_de_Y_con_ruido[,,2])
  arreglo_Y_gorro[,3,]<-t(arreglo_de_Y_con_ruido[,,3])
  arreglo_Y_gorro[,4,]<-t(arreglo_de_Y_con_ruido[,,4])
  #### transformando los errores
  ## se toma arreglo de transformaciones que sale de la función 
  max_transformacion<-sqrt(max(arreglo_de_transformaciones^2))
  arreglo_Y<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
  for (j in 1:n_len_t)
  {
    arreglo_Y[j,,]<-arreglo_de_transformaciones[,,j]%*%arreglo_Y_gorro[j,,]
  }
  ### dividiendo el arreglo entre la máxima entrada para que todas las entradas queden 
  ### con errores entre -1 y 1
  arreglo_Y<-arreglo_Y/max_transformacion
  datos<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
  datos[,1,]<-arreglo_de_fun_media_1 + arreglo_Y[,1,]
  datos[,2,]<-arreglo_de_fun_media_2 + arreglo_Y[,2,]
  datos[,3,]<-arreglo_de_fun_media_3 + arreglo_Y[,3,]
  datos[,4,]<-arreglo_de_fun_media_4 + arreglo_Y[,4,]
  
  #################################################################################################
  ## gráficos errores
  #################################################################################################
  # gráfico proceso 1 con dependencia
graficar_errores<-function()
{
  layout(matrix(c(1,2,3,4),nrow = 1))
  plot(t,arreglo_Y[,1,1],"l",col=1,main="",ylab="",ylim=c(-1,1))# errores proceso 1 con dependencia
  for (i in 1:N_tamano)
  {
    lines(t,arreglo_Y[,1,i],"l",col=gray(1-0.8*(i/N_tamano)))
  }
  plot(t,arreglo_Y[,2,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#errores proceso 2 con dependencia
  for (i in 1:N_tamano)
  {
    lines(t,arreglo_Y[,2,i],"l",col=gray(1-0.8*(i/N_tamano)))
  }
  plot(t,arreglo_Y[,3,1],"l",col=1,main="",ylab="",ylim=c(-1,1))#errores proceso 3 con dependencia
  for (i in 1:N_tamano)
  {
    lines(t,arreglo_Y[,3,i],"l",col=gray(1-0.8*(i/N_tamano)))
  }
  plot(t,arreglo_Y[,4,1],"l",col=1,main="",ylab="",ylim=c(-1,1))# errores proceso 1 con dependencia
  for (i in 1:N_tamano)
  {
    lines(t,arreglo_Y[,4,i],"l",col=gray(1-0.8*(i/N_tamano)))
  }    
}
#graficar_errores()
  #################################################################################################
  ## gráficos procesos 
  #################################################################################################
  # gráfico proceso 1 con dependencia
grafico_procesos_simulados<-function()
{
  y_max=max(datos[,1,])
  y_min=min(datos[,1,])
  plot(t,datos[,1,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X1",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,1,i],"l",col=rgb(escala_color,escala_color,1))
  }
  y_max=max(datos[,2,])
  y_min=min(datos[,2,])
  plot(t,datos[,2,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X2",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,2,i],"l",col=rgb(escala_color,escala_color,1))
  }
  y_max=max(datos[,3,])
  y_min=min(datos[,3,])
  plot(t,datos[,3,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X3",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,3,i],"l",col=rgb(escala_color,escala_color,1))
  }
  y_max=max(datos[,4,])
  y_min=min(datos[,4,])
  plot(t,datos[,4,1],"l",col=rgb(0.2,0.2,1),main="",ylab="X4",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,4,i],"l",col=rgb(escala_color,escala_color,1))
  }
}
#grafico_procesos_simulados()



#################################################################################################
#################################################################################################
#### TERMINA LA SIMULACIÓN DE LOS DATOS
#################################################################################################
#################################################################################################
#### EMPIEZA EL ANÁLISIS
#################################################################################################
#################################################################################################


#########################################################################################################
### encontrando componentes principales
lista_com_principales<-componentes_principales_multivariadas_continuas(datos)
lista_com_principales$array_mat_covarianzas
# valores propios
valores_propios<-lista_com_principales$array_valores_pro
graficar_val_propios_log<-function()
{
  layout(1)
  plot(t,log(valores_propios[1,]),"l",col=1,main="",ylab="",
       ylim=c(min(log(valores_propios)),max(log(valores_propios))))
  lines(t,log(valores_propios[2,]),col=2)
  lines(t,log(valores_propios[3,]),col=3)
  lines(t,log(valores_propios[4,]),col=4)
}
#graficar_val_propios_log()
# haciendo el cambio de base que hace las matrices de covarianza sean identidades
lista_transformados<-cambio_de_base_comp_principales(datos)
medias<-lista_transformados$arreglo_vec_medias

# gráficos media teórica y media estimada
graficar_medias<-function()
{
  layout(matrix(c(1,2,3,4),nrow = 1))
  y_min=min(medias[,1])
  y_max=max(medias[,1])
  plot(t,fun_media_1_antes,"l",col="green",main="Medias proceso 1",ylim=c(y_min,y_max))
  lines(t,fun_media_1_despues,"l",col="blue")
  lines(t,medias[,1],"l",col="red")
  y_min=min(medias[,2])
  y_max=max(medias[,2])
  plot(t,fun_media_2_antes,"l",col="green",main="Medias proceso 2",ylim=c(y_min,y_max))
  lines(t,fun_media_2_despues,"l",col="blue")
  lines(t,medias[,2],"l",col="red")
  y_min=min(medias[,3])
  y_max=max(medias[,3])
  plot(t,fun_media_3_antes,"l",col="green",main="Medias proceso 3",ylim=c(y_min,y_max))
  lines(t,fun_media_3_despues,"l",col="blue")
  lines(t,medias[,3],"l",col="red")
  y_min=min(medias[,4])
  y_max=max(medias[,4])
  plot(t,fun_media_4_antes,"l",col="green",main="Medias proceso 4",ylim=c(y_min,y_max))
  lines(t,fun_media_4_despues,"l",col="blue")
  lines(t,medias[,4],"l",col="red")
}
#graficar_medias()
# proceso transformado
proceso_transformado<-lista_transformados$nuevo_proceso
# las matrices de covarianza en cada punto deberían ser la identidad y los valores propios todos 
# deberían ser 1
prueba_comp<-componentes_principales_multivariadas_continuas(proceso_transformado)
prueba_comp$array_mat_covarianzas
prueba_comp$array_valores_pro
# gráficas  
graficar_transformados<-function()
{
  layout(matrix(c(1,2,3,4),nrow = 1))
  plot(t,proceso_transformado[,1,1],"l",main="",ylab="",
       ylim=c(min(proceso_transformado[,1,]),max(proceso_transformado[,1,])))
  for (i in 1:N_tamano)
  {
    lines(t,proceso_transformado[,1,i],col=i) 
  }  
  plot(t,proceso_transformado[,2,1],"l",main="",ylab="",
       ylim=c(min(proceso_transformado[,2,]),max(proceso_transformado[,2,])))
  for (i in 1:N_tamano)
  {
    lines(t,proceso_transformado[,2,i],col=i) 
  }  
  plot(t,proceso_transformado[,3,1],"l",main="",ylab="",
       ylim=c(min(proceso_transformado[,3,]),max(proceso_transformado[,3,])))
  for (i in 1:N_tamano)
  {
    lines(t,proceso_transformado[,3,i],col=i) 
  }  
  plot(t,proceso_transformado[,4,1],"l",main="",ylab="",
       ylim=c(min(proceso_transformado[,4,]),max(proceso_transformado[,4,])))
  for (i in 1:N_tamano)
  {
    lines(t,proceso_transformado[,4,i],col=i) 
  }
}
#graficar_transformados()
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
#
smooth_data_4<-Data2fd(argvals = t,y=proceso_transformado[,4,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_4,main="Datos suavizados")
#########################################################################################################
### análisis de punto de cambio
#########################################################################################################
# d no. de componentes principales
d1<-7
d2<-7
d3<-7
d4<-7
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

#
comp_principales_4<-pca.fd(smooth_data_4,nharm = d4,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_4$varprop
sum(comp_principales_4$varprop)
T_estadistico_4<-estadistico_T_N_x(comp_principales_4$scores,comp_principales_4$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3+T_estadistico_4
layout(1)
plot(seq(1,N_tamano,1)/N_tamano,T_estadistico,main="","l",ylab="",xlab="")
max(T_estadistico)
which.max(T_estadistico)




sum(T_estadistico)/N_tamano ### estadístico S
which.max(T_estadistico) # punto de cambio

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#### comparación sin hacer transformaciones. No recomendado
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
#
smooth_data_4<-Data2fd(argvals = t,y=datos[,4,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_4,main="Datos suavizados")
#########################################################################################################
### análisis de punto de cambio
#########################################################################################################
# d no. de componentes principales
d1<-7
d2<-7
d3<-7
d4<-7
comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_1$varprop
sum(comp_principales_1$varprop)
T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

#
comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2$varprop
sum(comp_principales_2$varprop)
T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
#
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3$varprop
sum(comp_principales_3$varprop)
T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)
#
comp_principales_4<-pca.fd(smooth_data_4,nharm = d4,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_4$varprop
sum(comp_principales_4$varprop)
T_estadistico_4<-estadistico_T_N_x(comp_principales_4$scores,comp_principales_4$values,FALSE)

#
T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3+T_estadistico_4
plot(seq(1,N_tamano,1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico)
which.max(T_estadistico)


#

sum(T_estadistico)/N_tamano
which.max(T_estadistico)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
### Haciendo deteccion de cambio  con MFPCA. No recomendado
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### Creating a multifunData object with 2 observations on the same domain
# Univariate elements

p1 <- funData(t,t(datos[,1,]))
plot(p1)
#
p2 <- funData(t,t(datos[,2,]))
plot(p2)
#
p3 <- funData(t,t(datos[,3,]))
plot(p3)
#
p4 <- funData(t,t(datos[,4,]))
plot(p4)

datos_multi_fun_data <- multiFunData(p1,p2,p3,p4)

uniExpansions <- list(list(type = "uFPCA", npc = 7), 
                      list(type = "uFPCA", npc = 7),
                      list(type = "uFPCA", npc = 7),
                      list(type = "uFPCA", npc = 7)) 

MFPCA_datos<- MFPCA(datos_multi_fun_data, M = 4, uniExpansions = uniExpansions)
summary(MFPCA_datos)
T_estadistico<-estadistico_T_N_x(MFPCA_datos$scores,MFPCA_datos$values,TRUE)
sum(T_estadistico)/N_tamano
which.max(T_estadistico)

