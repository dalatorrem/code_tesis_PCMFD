script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script.dir)
source("funciones_tesis_version_2.R")


datos<-read.csv("calidad_aire_con_faltantes.csv",sep=";")
names(datos)[1]<-"NO"
### faltantes
faltantes=0
for (j in 1:3)
{
  for (i in 1:(dim(datos)[1]))
  {
    if (is.na(datos[i,j]))
    {
      faltantes=faltantes+1
    }
  }
}
faltantes
for (j in 1:3)
{
  for (i in 1:(dim(datos)[1]))
  {
    if (is.na(datos[i,j]))
    {
      datos[i,j]<-(datos[(i-1),j]+datos[(i+1),j])/2 
    }
  }
}
for (j in 1:3)
{
  for (i in 1:(dim(datos)[1]))
  {
    if (is.na(datos[i,j]))
    {
      if (i>=168 & i<=(dim(datos)[1]-168))
      {
        if (!(is.na(datos[(i-168),j]) | is.na(datos[(i+168),j])))
        {
          datos[i,j]<-(datos[(i-168),j]+datos[(i+168),j])/2
        }
        else
        {
          if (is.na(datos[(i-168),j]))
          {
            datos[i,j]<-datos[(i+168),j]
          }
          if (is.na(datos[(i+168),j]))
          {
            datos[i,j]<-datos[(i-168),j]
          }
        }
        
      }
      else
      {
        if (i<168)
        {
          datos[i,j]<-datos[(i+168),j] 
        }
        if (i>(dim(datos)[1]-168)) 
        {
          datos[i,j]<-datos[(i-168),j]
        }
      }
      
    }
  }
}
for (j in 1:3)
{
  for (i in 1:(dim(datos)[1]))
  {
    if (is.na(datos[i,j]))
    {
      print(i)
    }
  }
}
datos<-datos[,1:2]  ## NOX es la suma de NO y NO2
write.table(datos,'calidad_aire_sin_faltantes.csv',sep = ';',row.names = FALSE)
View(datos)
no_tiempos<-24
t<-(1:no_tiempos)/no_tiempos
p=dim(datos)[2]
n=(dim(datos)[1])/no_tiempos
procesos<-array((1:(p*n*no_tiempos)),c(no_tiempos,p,n))
for (j in 1:p)
{
  procesos[,j,]<-matrix(datos[,j],nrow=n)
}
no_bases<-15

lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
procesos<-lista_cambio_de_base$nuevo_proceso

Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)

suavizamiento_1<-Data2fd(argvals = t,procesos[,1,],basisobj=Bfourier,dfscale = "gcv")
suavizamiento_2<-Data2fd(argvals = t,procesos[,2,],basisobj=Bfourier,dfscale = "gcv")

d1<-6
d2<-6

comp_principales_1<-pca.fd(suavizamiento_1,nharm = d1,
                           harmfdPar = fdPar(Bfourier),centerfns = TRUE)
comp_principales_2<-pca.fd(suavizamiento_2,nharm = d2,
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
(sum(T_estadistico)/length(T_estadistico))


