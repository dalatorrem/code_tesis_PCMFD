setwd("G:/Mi unidad/TESIS_MAESTRIA/programas R/definitivos/presentado")
source("funciones_tesis_version_2.R")
library(expm)
set.seed(1234)
y_1<-rnorm(5000,0,1)
set.seed(5678)
y_2<-rnorm(5000,0,1)
covarianza_teorica<-matrix(c(1,0.9,0.9,1),ncol=2)
raiz_cov_teorica<-sqrtm(covarianza_teorica)
# datos_y las componentes son independientes
datos_y<-matrix(c(y_1,y_2),ncol=5000,byrow = TRUE)
# agregando la dependencia
datos_x<-raiz_cov_teorica%*%datos_y
# covarianza estimada datos_y, el teórico es la identidad
cov(t(datos_y))
# covarianza estimada datos_x, el teórico es covarianza_teorica
cov(t(datos_x))
componentes_principales<-eigen(cov(t(datos_x)))
maximo<-sqrt(max(datos_x^2))
plot(datos_x[1,],datos_x[2,],xlim = c(-maximo-3,maximo+3),
     ylim=c(-maximo-3,maximo+3),"p",lwd=0.1,cex=0.6,
     xlab="X1",ylab="X2")
points(x=componentes_principales$vectors[1,2]*sqrt(componentes_principales$values[2])*6,
       y=componentes_principales$vectors[2,2]*sqrt(componentes_principales$values[2])*6
       ,col="green",lwd=2,cex=1)
points(x=componentes_principales$vectors[1,1]*sqrt(componentes_principales$values[1])*6,
       y=componentes_principales$vectors[2,1]*sqrt(componentes_principales$values[1])*6
       ,col="blue",lwd=2,cex=1)
# matriz para transformar los datos y eliminar la influencia de la covarianza
raiz_cov_mue_inv<-
  solve(componentes_principales$vectors%*%sqrt(diag(componentes_principales$values)))
datos_x_trans<-(raiz_cov_mue_inv%*%datos_x)
ext<-raiz_cov_mue_inv%*%componentes_principales$vectors
plot(datos_x_trans[1,],datos_x_trans[2,],xlim = c(-maximo-3,maximo+3),
     ylim=c(-maximo-3,maximo+3),"p",lwd=0.1,cex=0.6,
     xlab="Y1",ylab="Y2")
points(x=ext[1,2]*sqrt(componentes_principales$values[2])*6,
       y=ext[2,2]*sqrt(componentes_principales$values[2])*6
       ,col="green",lwd=2,cex=1)
points(x=ext[1,1]*sqrt(componentes_principales$values[1])*6,
       y=ext[2,1]*sqrt(componentes_principales$values[1])*6
       ,col="blue",lwd=2,cex=1)
cov(t(datos_x_trans))
######################################################################################
######################################################################################
######################################################################################


t_1<-seq(0,1,0.2)
t_2<-c(0.2,0.4,0.45,0.6,0.7,1)
set.seed(12)
error1<-rnorm(6,0,0.2)
f_1<-sin(2*pi*t_1)+error1
set.seed(1234578)
error2<-rnorm(6,0,0.2)
f_2<-sin(2*pi*t_2)+error2
plot(t_1,f_1,xlab="t",ylab="X",col="blue",lwd=4)
points(t_2,f_2,col="green",lwd=4)


#######################################################################################
#######################################################################################
# diagramas de dispersión en el tiempo

t<-seq(0,1,0.2)
p1<-sin(t*pi/2)
p2<-cos(t*pi/2)
p3<--cos(t*pi/2)*0.2
p4<-sin(t*pi/2)
no_N=20
mat_componentes<-matrix(c(p1,p2,p3,p4),nrow=6)
## verificar que los valores propios de A*A' no son cero y que las funciones de estos no se cruzan
arreglo_de_transformaciones<-verif_arreg_de_transformaciones(mat_componentes,t)
set.seed(12344)
normales<-rnorm(no_N*12,0,1)
arreglo_X_gorro<-array(normales,dim=c(6,2,no_N))
arreglo_X<-array(rep(0,no_N*12),dim=c(6,2,no_N))
for (j in 1:6)
{
  arreglo_X[j,,]<-arreglo_de_transformaciones[,,j]%*%arreglo_X_gorro[j,,]
}
Conf3x2 = matrix(c(1:6), nrow=2, byrow=TRUE)
layout(Conf3x2)
layout.show(6)
color=0.9
gris=color-color*(1:no_N)/no_N
for (j in 1:6) {
  tiempo=(j-1)/5
  titulo<-paste("t=",tiempo)
  plot(x=arreglo_X[j,1,],
       y=arreglo_X[j,2,],
       main = titulo,xlab=expression("X"[1]),ylab=expression("X"[2]),lwd=0.2,xlim=c(-5,5),
       ylim=c(-5,5),col=gray(gris),pch = 15)
}


Conf1x2 = matrix(c(1:2), nrow=1, byrow=TRUE)
layout(Conf1x2)
layout.show(2)
t_1=seq(0,1,0.01)
plot(x=t_1,y=predict(smooth.spline(t,arreglo_X[,1,1]),t_1)$y,col=gray(color-color/no_N),xlab="t",ylab=expression("X"[1]),"l",ylim=c(-5,5))
for (i in 2:no_N)
{
  lines(x=t_1,y=predict(smooth.spline(t,arreglo_X[,1,i]),t_1)$y,col=gray(color-i*color/no_N))
}
plot(x=t_1,y=predict(smooth.spline(t,arreglo_X[,2,1]),t_1)$y,col=gray(color-color/no_N),xlab="t",ylab=expression("X"[2]),"l",ylim=c(-5,5))
for (i in 2:no_N)
{
  lines(x=t_1,y=predict(smooth.spline(t,arreglo_X[,2,i]),t_1)$y,col=gray(color-i*color/no_N))
}


