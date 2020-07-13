direccion_script<-dirname(rstudioapi::getActiveDocumentContext()$path)
source("funciones_tesis_version_2.R")
t<- seq(0,1,0.01)



#####################################################################################################
##############
##### SE PRUEBA EL ENFOQUE SOBRE LA DERIVADA
set.seed(123456)
tamano_muestra_antes<-95
tamano_muestra_despues<-64

## usando simulación de procesos en Happ
Valtype<-"linear"
M_val_propios_base<-2
N_tamano<-tamano_muestra_antes+tamano_muestra_despues
objeto_de_Y_sin_ruido<- simFunData(N =N_tamano, 
                                   argvals = t,eFunType = "Fourier",eValType = Valtype, 
                                   M = M_val_propios_base)$simData


# Extrayendo el valor a mayor distancia del 0 de los simulados de Y 
max_dist_Y_sin_ruido<-sqrt(max((objeto_de_Y_sin_ruido@X)^2))
## fijando como desviacion el max extraido
desv<-max_dist_Y_sin_ruido
arreglo_de_Y_con_ruido<-t(addError(objeto_de_Y_sin_ruido,sd=desv)@X)





### se deja que el maximo error sea 1 solamente por convención

arreglo_de_Y_con_ruido<-arreglo_de_Y_con_ruido/(sqrt(max(arreglo_de_Y_con_ruido^2)))
boxplot(matrix(arreglo_de_Y_con_ruido,ncol = 1))
summary(matrix(arreglo_de_Y_con_ruido,ncol = 1))

## construyendo arreglo de medias de acuerdo al tamaño de la muestra 
####### Función antes del cambio
funcion_media_antes<-4*sin(2*pi*t)+8
arreglo_de_medias_antes<-kronecker(t(rep(1,tamano_muestra_antes)),funcion_media_antes)
####### Función después del cambio
### valor del cambio en el coseno 11
v_cambio<-0.4  #0.1   0.3  0.4  0.5 0.8
funcion_media_despues<-4*sin(2*pi*t)+8+v_cambio*cos(11*2*pi*t)    
arreglo_de_medias_despues<-kronecker(t(rep(1,tamano_muestra_despues)),funcion_media_despues)

plot(t,funcion_media_despues,main="",ylab="","l",col=2)
lines(t,funcion_media_antes,col=1)
## función de medias

arreglo_de_medias<-cbind(arreglo_de_medias_antes,arreglo_de_medias_despues)

### arreglo de simulaciones del proceso X 
arreglo_de_X<-arreglo_de_medias+arreglo_de_Y_con_ruido


# suavizando datos
no_bases<-35
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
plot(Bfourier,main="Gráfico elementos de la base")
smooth_data<-Data2fd(argvals = t,y=arreglo_de_X,basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data,main="Datos suavizados")

#Transformación de los datos, se modifica el elemento de la lista correspondiente a los coeficientes

smooth_data_derivada<-smooth_data
smooth_data_derivada$coefs<-matrix_derivada(no_bases)%*%as.matrix(smooth_data$coefs)
plot(smooth_data_derivada,main="Datos suavizados")




### detección de punto de cambio tradicional
d_no_comp_prin<-12
comp_principales<-pca.fd(smooth_data,nharm = d_no_comp_prin,harmfdPar = fdPar(Bfourier),centerfns = TRUE)

sum(comp_principales$varprop)
T_estadistico<-estadistico_T_N_x(comp_principales$scores,comp_principales$values,TRUE)
sum(T_estadistico)/N_tamano  # el mismo valor de estadistico S_N_D
estadistico_S_N_D(comp_principales$scores,comp_principales$values)
max(T_estadistico)
which.max(T_estadistico)


### detección de punto de cambio en la derivada
d_no_comp_prin<-15
comp_principales<-pca.fd(smooth_data_derivada,nharm = d_no_comp_prin,harmfdPar = fdPar(Bfourier),centerfns = TRUE)

sum(comp_principales$varprop)
T_estadistico<-estadistico_T_N_x(comp_principales$scores,comp_principales$values,TRUE)
sum(T_estadistico)/N_tamano  # el mismo valor de estadistico S_N_D
estadistico_S_N_D(comp_principales$scores,comp_principales$values)
max(T_estadistico)
which.max(T_estadistico)
