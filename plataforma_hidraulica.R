direccion_script<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(direccion_script)
source("funciones_tesis_version_2.R")
perfil <- read.table("profile.txt",header=FALSE)

filas<-perfil$V1==100 & perfil$V2==100 & perfil$V3==0 & perfil$V4==130
filas<-which(filas==TRUE)
filas
perfil[1665:1680,]





######## análisis univariado para Temperatura sensor 1

temp_s_1<-read.table("TS1.txt",header=FALSE)
temp_s_1<-temp_s_1[1665:1680,]
plot(x=1:60,temp_s_1[1,],type="l",col=1)
for (i in 1:dim(temp_s_1)[1]){
  i_tit=as.character(i)
  lines(x=1:60,temp_s_1[i,],type="l",main=i_tit,col=i)
}


no_bases=15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)

datos_temp_s_1<-t(temp_s_1)
t=0:59/60
temp_s_1_suavizada<-Data2fd(argvals = t,datos_temp_s_1,
                                         basisobj=Bfourier,dfscale = "gcv")

plot(temp_s_1_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylab=c("Temperatura °C"))

d_no_comp_prin<-6
comp_principales_temp_s_1<-pca.fd(temp_s_1_suavizada,
                                             nharm = d_no_comp_prin,
                                             harmfdPar = fdPar(Bfourier),
                                             centerfns = TRUE)
comp_principales_temp_s_1$varprop
sum(comp_principales_temp_s_1$varprop)
est_s<-estadistico_S_N_D(comp_principales_temp_s_1$scores,
                         comp_principales_temp_s_1$values)

estT<-estadistico_T_N_x(comp_principales_temp_s_1$scores,
                        comp_principales_temp_s_1$values,TRUE)

max(estT)
which.max(estT)
sum(estT)/16
est_s


######## análisis univariado para Temperatura sensor 2

temp_s_2<-read.table("TS2.txt",header=FALSE)
temp_s_2<-temp_s_2[1665:1680,]
plot(x=1:60,temp_s_2[1,],type="l",col=1)
for (i in 1:dim(temp_s_2)[1]){
  i_tit=as.character(i)
  lines(x=1:60,temp_s_2[i,],type="l",main=i_tit,col=i)
}


no_bases=15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)

datos_temp_s_2<-t(temp_s_2)
t=0:59/60
temp_s_2_suavizada<-Data2fd(argvals = t,datos_temp_s_2,
                            basisobj=Bfourier,dfscale = "gcv")

plot(temp_s_2_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylab=c("Temperatura °C"))

d_no_comp_prin<-6
comp_principales_temp_s_2<-pca.fd(temp_s_2_suavizada,
                                 nharm = d_no_comp_prin,
                                 harmfdPar = fdPar(Bfourier),
                                 centerfns = TRUE)
comp_principales_temp_s_2$varprop
est_s<-estadistico_S_N_D(comp_principales_temp_s_2$scores,
                         comp_principales_temp_s_2$values)

estT<-estadistico_T_N_x(comp_principales_temp_s_2$scores,
                        comp_principales_temp_s_2$values,TRUE)

max(estT)
which.max(estT)
est_s
sum(estT)/length(estT)

######## análisis univariado para Temperatura sensor 3

temp_s_3<-read.table("TS3.txt",header=FALSE)
temp_s_3<-temp_s_3[1665:1680,]
plot(x=1:60,temp_s_3[1,],type="l",col=1)
for (i in 1:dim(temp_s_3)[1]){
  i_tit=as.character(i)
  lines(x=1:60,temp_s_3[i,],type="l",main=i_tit,col=i)
}


no_bases=15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)

datos_temp_s_3<-t(temp_s_3)
t=0:59/60
temp_s_3_suavizada<-Data2fd(argvals = t,datos_temp_s_3,
                            basisobj=Bfourier,dfscale = "gcv")

plot(temp_s_3_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylab=c("Temperatura °C"))

d_no_comp_prin<-6
comp_principales_temp_s_3<-pca.fd(temp_s_3_suavizada,
                                  nharm = d_no_comp_prin,
                                  harmfdPar = fdPar(Bfourier),
                                  centerfns = TRUE)
comp_principales_temp_s_3$varprop
est_s<-estadistico_S_N_D(comp_principales_temp_s_3$scores,
                         comp_principales_temp_s_3$values)

estT<-estadistico_T_N_x(comp_principales_temp_s_3$scores,
                        comp_principales_temp_s_3$values,TRUE)

max(estT)
which.max(estT)
sum(estT)/16
est_s

######## análisis univariado para Temperatura sensor 4

temp_s_4<-read.table("TS4.txt",header=FALSE)
temp_s_4<-temp_s_4[1665:1680,]
plot(x=1:60,temp_s_4[1,],type="l",col=1)
for (i in 1:dim(temp_s_4)[1]){
  i_tit=as.character(i)
  lines(x=1:60,temp_s_4[i,],type="l",main=i_tit,col=i)
}


no_bases=15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)

datos_temp_s_4<-t(temp_s_4)
t=0:59/60
temp_s_4_suavizada<-Data2fd(argvals = t,datos_temp_s_4,
                            basisobj=Bfourier,dfscale = "gcv")

plot(temp_s_4_suavizada,main="Temperaturas exteriores suavizadas",
     xlab="Transcurso del día",ylab=c("Temperatura °C"))

d_no_comp_prin<-5
comp_principales_temp_s_4<-pca.fd(temp_s_4_suavizada,
                                  nharm = d_no_comp_prin,
                                  harmfdPar = fdPar(Bfourier),
                                  centerfns = TRUE)
comp_principales_temp_s_4$varprop
sum(comp_principales_temp_s_4$varprop)
est_s<-estadistico_S_N_D(comp_principales_temp_s_4$scores,
                         comp_principales_temp_s_4$values)

estT<-estadistico_T_N_x(comp_principales_temp_s_4$scores,
                        comp_principales_temp_s_4$values,TRUE)

max(estT)
which.max(estT)
sum(estT)/16
est_s

#############################################################################################
###### metodología multivariada
#############################################################################################
dim(temp_s_1)
n<-dim(temp_s_1)[1]
no_tiempos<-dim(temp_s_1)[2]
dim(datos_temp_s_1) #verificando dimensiones :)
dim(datos_temp_s_2)      #verificando dimensiones :)
dim(datos_temp_s_3)      #verificando dimensiones :)
dim(datos_temp_s_4)      #verificando dimensiones :)
# generando arreglo 
# proceso1 temperatura exterior
# proceso2 temperatura hab
p=4
procesos<-array((1:(p*n*no_tiempos)),c(no_tiempos,p,n))
for (j in 1:n)
{
  procesos[,1,j]<-datos_temp_s_1[,j]
  procesos[,2,j]<-datos_temp_s_2[,j]
  procesos[,3,j]<-datos_temp_s_3[,j]
  procesos[,4,j]<-datos_temp_s_4[,j]
}



lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)




raiz_val_propios<-sqrt(lista_cambio_de_base$array_valores_pro)  

# raiz de los valores propios en el tiempo
ymin<-min(raiz_val_propios)
ymax<-max(raiz_val_propios)
plot(t,raiz_val_propios[1,],"l",col=1,ylim = c(ymin,ymax))
lines(t,raiz_val_propios[2,],"l",col=2)
lines(t,raiz_val_propios[3,],"l",col=3)
lines(t,raiz_val_propios[4,],"l",col=4)
proceso_transformado<-lista_cambio_de_base$nuevo_proceso




no_bases<-15
Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
#
smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_1,main="Componente 1 suavizamientos")
#
smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_2,main="Componente 2 suavizamientos")

smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_3,main="Componente 3 suavizamientos")

smooth_data_4<-Data2fd(argvals = t,y=proceso_transformado[,4,],basisobj=Bfourier,dfscale = "gcv")
plot(smooth_data_4,main="Componente 3 suavizamientos")


d1<-3
d2<-7
d3<-8
d4<-7

comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

comp_principales_4<-pca.fd(smooth_data_4,nharm = d4,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)

comp_principales_1$varprop
sum(comp_principales_1$varprop)

comp_principales_2$varprop
sum(comp_principales_2$varprop)


comp_principales_3$varprop
sum(comp_principales_3$varprop)

comp_principales_4$varprop
sum(comp_principales_4$varprop)




T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)

T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)

T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)

T_estadistico_4<-estadistico_T_N_x(comp_principales_4$scores,comp_principales_4$values,FALSE)
#



T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3+T_estadistico_4
plot(seq(1,length(T_estadistico[,1]),1),T_estadistico,main="Gráfico estadístico T","l")
max(T_estadistico) #11.81754
which.max(T_estadistico) #10
sum(T_estadistico/n) #6.273834




