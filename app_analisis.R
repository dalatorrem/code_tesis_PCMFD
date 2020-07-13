library(shiny)
library(shinythemes)
library(dplyr)
#### paquetes para datos funcionales
library(fda)
library(funData)
library(MFPCA)
library(readxl)
library(abind)
library(rmutil)
library(dplyr)
#### paquetes para datos funcionales

#################################################################################################
##funciones app



#########################################################################################################
#########################################################################################################
### funciones analisis
#########################################################################################################
#########################################################################################################
plotear_var_discretas<-function(datos,no_tiempos,no_bases,comp_prin)
{
  no_tiempos<-eval(parse(text=no_tiempos))
  datos<-traer_datos(datos,no_tiempos)
  procesos<-datos$procesos
  nom_var<-datos$nombres
  p=dim(procesos)[2]
  n=dim(procesos)[3]
  # preparar para graficos
  Conf= matrix(1:p, nrow=1)
  graf<-layout(Conf)
  layout.show(graf)
  t=(1:no_tiempos)
  # plots y arreglos
  y_min=min(procesos)
  y_max=max(procesos)
  for (j in 1:p)
  {
    y_min=min(procesos[,j,])
    y_max=max(procesos[,j,])
    plot(x=t,y=procesos[,j,1],col=gray(0.2),xlab="t",
         ylab=nom_var[j],"l",ylim=c(y_min,y_max))
    for (i in 2:n)
    {
      gris=0.9-0.7*(n-i)/n
      lines(x=t,y=procesos[,j,i],col=gray(gris))
    }
  }    
}
################################
traer_datos<-function(datos,no_tiempos)
{
  datos<-read.csv(datos$datapath,sep=";")
  p=dim(datos)[2]
  n=(dim(datos)[1])/no_tiempos
  procesos<-array((1:(p*n*no_tiempos)),c(no_tiempos,p,n))
  for (j in 1:p)
  {
    procesos[,j,]<-matrix(datos[,j],nrow=n)
  }
  nombres=names(datos)
  datos<-list(procesos=procesos,nombres=nombres)
  return(datos)
}
###############################################################################################################################################
# detectar
###############################################################################################################################################
detectar<-function(datos,no_tiempos,no_bases,com_prin)
{
  no_tiempos<-eval(parse(text=no_tiempos))
  com_prin<-eval(parse(text = com_prin))
  no_bases<-eval(parse(text = no_bases))
  datos<-traer_datos(datos,no_tiempos)
  procesos<-datos$procesos
  nom_var<-datos$nombres
  p=dim(procesos)[2]
  t=(1:no_tiempos)/no_tiempos
  n=dim(procesos)[3]
  ## aplicando componentes principales multivariadas
  Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
  if (p==1)
  {
    d=com_prin[1]
    suavizamiento<-Data2fd(argvals = t,procesos[,1,],basisobj=Bfourier,dfscale = "gcv")
    comp_principales<-pca.fd(suavizamiento,nharm = d,
                             harmfdPar = fdPar(Bfourier),centerfns = TRUE)
    secu_graf=c(1,2,3,4,4,4)
    Conf= matrix(secu_graf, nrow=2,byrow = TRUE)
    graf<-layout(Conf)
    layout.show(graf)
    barplot(cumsum(100*comp_principales$varprop), main="Variabilidad Explicada",xlab="No. de componentes",
            ylab="% var explicada",names.arg=1:d,ylim=c(0,101))
    abline( h =95,col="red")
    abline( h =99,col="red")
    ###### trabajando....
    T_estadistico<-estadistico_T_N_x(comp_principales$scores,comp_principales$values,FALSE)
  }
  else
  {
    # matrix para graficos
    secu_graf=c(1:p,rep(p+1,(p)),p+2,rep(p+3,(p-1)))
    Conf= matrix(secu_graf, nrow=3,byrow = TRUE)
    graf<-layout(Conf)
    layout.show(graf)
    #
    lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
    proceso_transformado<-lista_cambio_de_base$nuevo_proceso
    T_estadistico=0
    for (j in 1:p)
    {
      d=com_prin[j]
      suavizamiento<-Data2fd(argvals = t,y=proceso_transformado[,j,],basisobj=Bfourier,dfscale = "gcv")
      comp_principales<-pca.fd(suavizamiento,nharm = d,
                               harmfdPar = fdPar(Bfourier),centerfns = TRUE)
      barplot(cumsum(100*comp_principales$varprop), main="Variabilidad Explicada",xlab="No. de componentes",
              ylab="% var explicada",names.arg=1:d,ylim=c(0,101))
      abline( h =95,col="red")
      abline( h =99,col="red")
      ###### trabajando....
      T_estadistico<-T_estadistico+estadistico_T_N_x(comp_principales$scores,comp_principales$values,FALSE)  
    }
    d=sum(com_prin[1:p])
    ### grafico para ver los valores propios
    lista_cambio_de_base<-cambio_de_base_comp_principales(procesos)
    raiz_val_propios<-sqrt(lista_cambio_de_base$array_valores_pro)
    ymin<-min(raiz_val_propios)
    ymax<-max(raiz_val_propios)
    for (j in 1:p)
    {
      if (j==1)
      {
        plot(t,raiz_val_propios[1,],"p",col=1,
             main="Raiz los valores propios de B(t)B(t)'(estimados y suavizamiento)",ylim=c(ymin,ymax),ylab="")
        ajuste_val<-lm(raiz_val_propios[j,]~bs(t,round(no_tiempos/2)))
        lines(t,ajuste_val$fitted.values,col=1)
      }
      else
      {
        lines(t,raiz_val_propios[j,],"p",col=j)
        ajuste_val<-lm(raiz_val_propios[j,]~bs(t,round(no_tiempos/2)))
        lines(t,ajuste_val$fitted.values,col=j)
      }
      
    }
  }# cierra el else de mutivariado
  
  plot(seq(1,n,1),T_estadistico,
       main="Estadistico T
       (Posible punto de cambio)",
       ylab="","l",xlab="",ylim=c(0,1.7*max(T_estadistico)))
  points(which.max(T_estadistico),max(T_estadistico),pch=23,col="green",bg="red",cex = 3.5)
  coor=paste0("Candidato a PC en:",toString(which.max(T_estadistico)))
  text(x=n/2,y=1.3*max(T_estadistico),
       label= coor,cex = 1.4)
  estadistico_prueba=sum(T_estadistico)/length(T_estadistico)
  y_max=1.2*max(tabla_distri_kd,estadistico_prueba)
  x_max=max(tabla_distri_kd[,d],estadistico_prueba)
  plot(x=estadistico_prueba,y=d,pch=17,col="green",cex = 4,
       main="Distribucion Kd y Estadistico de prueba",
       ylab="Suma componentes",xlab="valor",
       ylim=c(d-1,d+4),xlim=c(0,x_max))
  abline(h =d,col=gray(0.7))
  points(x=tabla_distri_kd[2,d],y=d,pch=4,cex=4,col=gray(0.1))
  points(x=tabla_distri_kd[3,d],y=d,pch=8,cex=4,col=gray(0.1))
  points(x=tabla_distri_kd[1,d],y=d,pch=3,cex=4,col=gray(0.1),"p")
  
  legend( 0,d+4, c("Est. Prueba",NA,"Percentil 90", "Percentil 95","Percentil 99"),
          pch = c(17 ,NA,3, 4,8),col=c("green",NA,col=gray(0.1),col=gray(0.1),col=gray(0.1)))
  
}

##############################################################################################################################



### subindices:
# l      para no. de componentes principales    d  d_no_comp_prin
# i      para tamaño de la muestra              N
# j      para la longitud de t                  n  n_lent_t
# k      para dimensiones                       p
# n_b    para no. de bases                      no_bases
# sub_i  para 2^p                               2^p
# sub_cp para p^2                               p^2
#########################################################################################################
#########################################################################################################
## Estadistico S_d ver Berkes et al Detecting changes in the mean of functional observations pag 932, 3.4
#########################################################################################################
#########################################################################################################

# ingresa la matriz de scores despues de haber realizado componentes principales para datos 
# funcionales, el obejto que entra tiene la estructura de PCA$scores, se debe tener cuidado 
# con la estructura si se usa otra funcion de componentes principales funcionales

estadistico_S_N_D=function(scores,valores_propios)
{
  N<-dim(scores)[1] # tamaño de la muestra
  d_no_comp_prin<-dim(scores)[2]
  # prueba funcion 1
  s=0
  for (l in 1:d_no_comp_prin)
  {
    inv_val<-1/valores_propios[l]
    score_l<-scores[,l]
    suma=0
    for (i in 1:N)
    {
      sub_1<-score_l[1:i]
      a<-sum(sub_1) # primera suma dentro del parentesis que esta al cuadrado
      b<-sum(score_l) # segunda suma dentro del parentesis que esta al cuadrado
      parentesis_al_cuadrado<-(a-(i/N)*b)^2
      suma=suma+parentesis_al_cuadrado
    }
    s=s+inv_val*suma
  }
  s=s/N^2
  return(s)
}



#########################################################################################################
#########################################################################################################
## Estadistico T_d ver Berkes et al Detecting changes in the mean of functional observations pag 931, 3.3
#########################################################################################################
#########################################################################################################

# ingresa la matriz de scores despues de haber realizado componentes principales para datos 
# funcionales, el obejto que entra tiene la estructura de PCA$scores, se debe tener cuidado 
# con la estructura si se usa otra funcion de componentes principales funcionales

estadistico_T_N_x=function(scores,values,plot_bool)
{
  N<-dim(scores)[1]
  d_no_comp_prin<-dim(scores)[2]
  est_t_x<-c()
  for (i in 1:N)
  {
    suma=0
    for (l in 1:d_no_comp_prin) 
    {
      inv_val<-1/values[l]
      score_l<-scores[,l]
      sub_1<-score_l[1:i]
      a<-sum(sub_1) # primera suma dentro del parentesis que esta al cuadrado
      b<-sum(score_l) # segunda suma dentro del parentesis que esta al cuadrado
      parentesis_al_cuadrado<-(a-(i/N)*b)^2
      suma=suma+inv_val*parentesis_al_cuadrado
    }
    T_i<-suma/N
    est_t_x<-rbind(est_t_x,T_i)
  }
  if (plot_bool==TRUE)
  {
    plot(seq(1,N,1),est_t_x,main="Grafico estadistico T","l")
  }
  return(est_t_x)
}   



#########################################################################################################
#########################################################################################################
##  SE INGRESA UN VALOR K ENTERO Y DEVUELVE LA MATRIZ DE DERIVADA PARA UNA BASE DE FOURIER DE 
##  ORDEN K, LA BASE DE FOURIER ESTa EN EL ORDEN DADO POR LA FUNCIoN create.fourier.basis()$names
#########################################################################################################
#########################################################################################################



matrix_derivada<-function(no_bases)
{
  ceros<-t(as.matrix(rep(0,no_bases)))
  mat_derivada<-ceros
  for (n_b in 2:no_bases) 
  {
    # para las posiciones pares o sea para los senos lo derivada es la contraccion del seno (2pi)
    # y el signo es positivo, para el coseno es analogo salvo que el signo es negativo. para 1 der=0
    der<-floor(n_b/2)*2*pi*(-1)^n_b  
    derivada<-ceros
    # para las posiciones pares (senos) lo envia uno adelante o sea al coseno del mismo argumento
    # para las posiciones impares (cosenos) lo envia uno atras o sea al seno del mismo argumento
    posicion<-n_b+1-2*(n_b%%2) 
    derivada[1,posicion]<-der
    mat_derivada<-rbind(mat_derivada,derivada)
  } 
  return(mat_derivada)
}



#########################################################################################################
#########################################################################################################
#### Componentes principales multivariados
#########################################################################################################
#########################################################################################################

# el argumento es un proceso p-variado sin datos faltantes. Los intervalos en los que se tomaron
# las medidas tienen que ajustarse al intervalo [0,1] haciendo una traslacion y un reescalamiento
# del tiempo. 
# en esta funcion se calculan los componentes principales p-variados en cada division de tiempo
# del interavlo [0,1]
# el argumento procesos es un arreglo [len,p,N]
# len<- numero de intervalos en [0,1]
# p<- dimensiones del proceso (p-variado)
# N<- tamaño de la muestra
# la funcion retorna una lista en la que hay:
# arreglos de matrices de covarianza pxp
# matrices con los p-valores propios en cada tiempo
# matriz ortonormal con los vectores propios buscadas para que sean continuas

componentes_principales_multivariadas_continuas<-function(procesos)
{
  n_len_t<-length(procesos[,1,1])
  p <-length(procesos[1,,1])
  N<-length(procesos[1,1,])
  
  arreglo_mat_varianzas<-array(rep(0,(p*p*n_len_t)),c(p,p,n_len_t)) # creando un arreglo pxpxn_len_t para luego modificarlo
  arreglo_mat_vectores_propios<-array(rep(0,(p*p*n_len_t)),c(p,p,n_len_t)) # creando un arreglo pxpxn_len_t para luego modificarlo
  arreglo_valores_propios<-matrix(rep(0,(p*n_len_t)),ncol=n_len_t) # creando una matrix pxlen_t para luego modificarla
  
  for (j in 1:n_len_t)
  {
    varianza<-var(t(procesos[j,,]))# estimador de la varianza
    arreglo_mat_varianzas[,,j]<-varianza
    propios<-eigen(varianza)
    arreglo_valores_propios[,j]<-t(propios$values)
    if(j==1)
    {
      vectores_propios_t_ant=diag(p) # definiendo una Base ortonormal para seleccionar
      # los vectores propios lo mas cercano a esa base, dado que hay 2^p posibilidades para seleccionar 
      # la base ortogonal de vectores propios, en el paso 1 se selecciona como los vectores canonicos
    }else
    {
      vectores_propios_t_ant<-vectores_propios_t # para evitar al maximo cambios bruscos en 
      # la matriz ortonormal con la base de vectores propios se selecciona la base de vectores propios 
      # que este lo mas cerca posible a la anterior base 
      # para esto se define vectores_propios_t_ant
    }
    vectores_propios_t<-seleccion_vectores_propios_continuidad(vectores_propios_t_ant,propios$vectors)
    arreglo_mat_vectores_propios[,,j]<-vectores_propios_t
  }
  lista_retorno<-list(array_mat_vec_pro=arreglo_mat_vectores_propios,
                      array_mat_covarianzas=arreglo_mat_varianzas,
                      array_valores_pro=arreglo_valores_propios)
  return(lista_retorno)
}




#########################################################################################################
#########################################################################################################
#### Seleecion de la base de vectores propios para garantizar continuidad
#########################################################################################################
#########################################################################################################

# esta funcion recibe dos matrices que serian la base ortonormal en el instante de tiempo anterior y la 
# de t. La funcion hace todas las posibles de los signos de las columnas y devuelve aquella matrix que 
# tenga la menor distancia de frobenius a la base ortonormal en el instante anterior.


seleccion_vectores_propios_continuidad<-function(vectores_propios_t_ant,vectores_propios_t_des)
{
  p<-dim(vectores_propios_t_des)[1]  # dimension
  distancias<-c()                    #
  arreglo_matrices_vectores_propios<-array(rep(0,p*p*(2^p)),c(p,p,2^p))
  for (sub_i in 1:(2^p))
  {
    matriz_vectores_propios<-vectores_propios_t_des
    signos<-intToBits(sub_i-1)[1:p] # sistema binario mostrando las potencias desde 2^0 hasta 2^p del valor sub_i-1 
    for (k in 1:p)
    {
      
      caracteristica<-if(signos[k]==intToBits(0)[1]){0}else{1}
      matriz_vectores_propios[,k]<-((-1)^caracteristica)*matriz_vectores_propios[,k]
    }
    arreglo_matrices_vectores_propios[,,sub_i]<-matriz_vectores_propios
    # calculando la diferencia de las matrices
    cambio_vec_propios<-(matriz_vectores_propios-vectores_propios_t_ant)
    if(det(matriz_vectores_propios)<0)  # solo se tienen en cuenta matrices con determinante positivo, grupo lineal
    {
      distancia<-Inf
    }
    else
    {
      distancia<-norm(cambio_vec_propios,"F")
    }
    distancias<-c(distancias,distancia)
  }
  #print(arreglo_matrices_vectores_propios)
  posicion_min_distancia<-which.min(distancias)
  return(arreglo_matrices_vectores_propios[,,posicion_min_distancia])
}



#########################################################################################################
#########################################################################################################
#### Cambio de base basado en componentes principales
#########################################################################################################
#########################################################################################################

# Esta funcion recibe los procesos como un arreglo de la misma forma que las demas funciones
# En la funcion se realiza en cada instante de tiempo una transformacion que deje la matriz 
# de covarianza como la identidad
# returna un arreglo en el que hay:
# una matriz de medias 
# un arreglo con los procesos transformados

cambio_de_base_comp_principales<-function(procesos)
{
  lista_retorno<-componentes_principales_multivariadas_continuas(procesos)
  array_mat_vec_prop<-lista_retorno$array_mat_vec_pro
  array_valores_pro<-lista_retorno$array_valores_pro
  
  n_len_t<-length(procesos[,1,1])
  N<-length(procesos[1,1,])
  p<-length(procesos[1,,1])
  
  # definiendo arreglos con ceros
  arreglo_vec_medias<-matrix(rep(0,(p*n_len_t)),ncol=p)
  procesos_con_errores_correlacionados<-array(rep(0,(p*n_len_t*N)),c(n_len_t,p,N))
  procesos_con_errores_no_correlacionados_estandarizados<-array(rep(0,(p*n_len_t*N)),c(n_len_t,p,N))
  nuevo_proceso<-array(rep(0,(p*n_len_t*N)),c(n_len_t,p,N))
  
  
  for (j in 1:n_len_t)
  {
    # calculando medias p-variadas de las n-muestras en el tiempo indexado por j
    for (k in 1:p)
    {
      arreglo_vec_medias[j,k]<-mean((procesos[j,k,]))
    }
    
    #definiendo la matriz Q^-0.5, para transformar los datos, obviamente esta es la misma matriz para 
    # cada uno de las muestras en el tiempo t
    elementos_diagonal<-1/sqrt(array_valores_pro[,j])
    #print(elementos_diagonal)
    lambda_inv_raiz<-diag(elementos_diagonal)
    #print(lambda_inv_raiz)
    for (i in 1:N) 
    {
      errores_correlacionados<-procesos[j,,i]-arreglo_vec_medias[j,]
      procesos_con_errores_correlacionados[j,,i]<-errores_correlacionados
      
      #### Sea S una matriz de covarianza definida positiva, luego cumple S=UD^2U' donde D^2 es una 
      #### diagonal con los valores propios (distintos) y positivos  de S, y U es una matriz ortonormal con det=1.
      #### observe que una matriz A= UD^-1 genera la misma matriz de covarianza S
      #### por lo tanto Q^-0.5*V'=A*^-1 quita las correlaciones
      
      errores_no_correlacionados<-t(lambda_inv_raiz%*%t(array_mat_vec_prop[,,j])%*%errores_correlacionados)
      procesos_con_errores_no_correlacionados_estandarizados[j,,i]<-errores_no_correlacionados
      nuevo_proceso[j,,i]<-arreglo_vec_medias[j,]+errores_no_correlacionados
    }
  }
  lista_retorno[["arreglo_vec_medias"]]<-arreglo_vec_medias
  lista_retorno[["procesos_con_errores_correlacionados"]]<-procesos_con_errores_correlacionados
  lista_retorno[["procesos_con_errores_no_correlacionados_estandarizados"]]<-procesos_con_errores_no_correlacionados_estandarizados
  lista_retorno[["nuevo_proceso"]]<-nuevo_proceso
  return(lista_retorno)
}



#########################################################################################################
#########################################################################################################
#### funcion para verificar si las componentes de una matriz cumplen las condiciones de Berrendero  
#########################################################################################################
#########################################################################################################

verif_arreg_de_transformaciones<-function(matriz,t)
{
  p<-sqrt(dim(matriz)[2])
  n_len_t<-dim(matriz)[1]
  arreglo_de_transformaciones<-array(rep(0,p*p*n_len_t),c(p,p,n_len_t))
  for (j in 1:n_len_t)
  {
    mat_t<-c()
    for (sub_cp in 1:p^2)
    {
      mat_t<-c(mat_t,matriz[j,sub_cp])
    }
    matriz_t<-matrix(mat_t,ncol = p)
    arreglo_de_transformaciones[,,j]<-matriz_t
  }
  
  # arreglo de matrices de covarianza teoricos 
  
  arreglo_de_covarianzas<-array(rep(0,p*p*n_len_t),c(p,p,n_len_t))
  
  for (j in 1:n_len_t)
  {
    arreglo_de_covarianzas[,,j]<- arreglo_de_transformaciones[,,j]%*%t(arreglo_de_transformaciones[,,j])
  }
  
  # valores propios de las matrices de covarianza en cada instante del tiempo, estos valores propios 
  # deben cumplir las hipotesis de Berrendero, o sea los lambda en el tiempo deben ser continuos, 
  # deben ser una sucesion decreciente de funciones
  
  valores_propios<-matrix(rep(0,p*n_len_t),ncol = n_len_t)
  for (j in 1:n_len_t)
  {
    valores_propios[,j]<-eigen(arreglo_de_covarianzas[,,j])$values
  }
  
  
  # graficos de las  funciones de los valores propios
  
  
  #ymax<-max((valores_propios))
  #ymin<-min((valores_propios))
  
  #plot(t,(valores_propios[1,]),ylim = c(ymin,ymax),"l",col=1, 
  #     main="",ylab="")# grafico de valores propios
  #for (k in 2:p)
  #{
  #  lines(t,valores_propios[k,],col=k)
  #}
  
  
  # graficos de los logaritmos de las  funciones de los valores propios
  
  
  ymax<-max(log(valores_propios))
  ymin<-min(log(valores_propios))
  
  plot(t,log(valores_propios[1,]),ylim = c(ymin,ymax),"l",col=1,
       main="Logaritmo natural  de los valores propios de B(t)B(t)' ",ylab="") # Grafico de logaritmos naturales de los valores propios en el intervalo
  for (k in 2:p)
  {
    lines(t,log(valores_propios[k,]),col=k)
  }
  return(arreglo_de_transformaciones)
}
tabla_distri_kd<-matrix(c(
  0.345165,0.606783,0.842567,1.065349,1.279713,1.4852,1.690773,1.897365,2.096615,2.288572,2.496635,2.686238,2.884214,3.066906,3.268958,3.462039,3.650724,3.837678,4.024313,4.2148,4.404677,4.591972,4.778715,4.965613,5.159057,5.346543,5.521107,5.714145,5.885108,6.083306,
  0.460496,0.748785,1.00139,1.239675,1.469008,1.684729,1.895557,2.124153,2.322674,2.526781,2.744438,2.949004,3.147604,3.336262,3.544633,3.740248,3.949054,4.136169,4.327286,4.532917,4.718904,4.908332,5.101896,5.303462,5.495721,5.688849,5.866095,6.068351,6.24277,6.444772,
  0.740138,1.072101,1.352099,1.626695,1.866702,2.12595,2.342252,2.589244,2.809778,3.033944,3.268031,3.491102,3.708033,3.903995,4.116829,4.317087,4.55465,4.734714,4.974172,5.156282,5.369309,5.576596,5.759427,5.973941,6.203718,6.393582,6.572949,6.771058,6.977607,7.186491
),nrow = 3,byrow = TRUE)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  #######################################################################
  # titulo
  #######################################################################
  titlePanel("ANALISIS: DETECCION DE PUNTO
             DE CAMBIO EN LA MEDIA PARA DATOS FUNCIONALES MULTIVARIADOS"),
  sidebarLayout(
    sidebarPanel(style = "overflow-y:scroll;max-height: 90vh;  position:relative;",
                 #######################################################################
                 # Procedimientos
                 #######################################################################
                 h6(strong("Cargar datos:")),
                 helpText(h6("Se debe cargar un archivo csv con separador punto y coma ';', cada variable debe ser una columna (max 5) 
                          con su respectivo nombre de columna, y cada fila 
                          corresponde a un tiempo en el que se realizo la medicion. si por ejemplo, los datos corresponden a 
                          mediciones de temperaturas tomadas cada hora durante 365 dias y cada dia sera analizado como un dato funcional,
                          se tendra n=365, se debe especificar el valor d=24, y el numero de filas del csv deberia ser igual a 365*24. 
                          no se admiten valores faltantes. introduzca el valor de d y luego cargue los datos.")),
                 textInput(inputId = "d",label = "",value = 60),
                 fileInput(inputId = "datos",label =  "",accept =".csv",buttonLabel = "Cargar...",
                           placeholder = "ej: datos_hidra_app.csv"),
                 h6(strong("Numero de Funciones de la base:")),
                 helpText(h6("En esta aplicacion se esta trabajando con bases de fourier.
                              el parametro de regularizacion es obtenido por GCV.")),
                 textInput(inputId="no_bases", label="", value = 25),
                 h6(strong("Numero de componetes:")),
                 helpText(h6("Coloque un arreglo con el numero de componentes principales para el acp univariado funcional de cada variable.
                            cada valor debe ser menor que el numero de bases.  Es recomendable escoger  un numero de componentes que logre
                            explicar entre el 95% y el 99%.  En esta aplicacion se puede decidir si la suma de 
                            los  valores es menor o igual a 30. 
                            Tenga en cuenta el grafico: estimado de la variabilidad explicada")),
                 textInput(inputId="comp_prin", label="", value ="c(5,5,5,5,5)")
                 
                 
                 
    ),
    mainPanel(style = "overflow-y:scroll;max-height: 90vh;  position:relative;",
              fluidRow(
                h6(strong("Observaciones Discretizadas (sin suavizamientos)")),
                plotOutput(outputId="plot_var_discretas",height = "200px"),
                plotOutput(outputId="plot_analisis",height = "500px")
              )
    )
  ))

# Define server ----
server <- function(input, output) {
  output$plot_var_discretas <- renderPlot({
    if (is.null(input$datos))
      {return(NULL)}
    else
      {plotear_var_discretas(input$datos,input$d)}
  })
  output$plot_analisis <- renderPlot({
    if (is.null(input$datos))
    {return(NULL)}
    else
    {detectar(input$datos,input$d,input$no_bases,input$comp_prin)}
  })
}


# Create Shiny app and deploy ----
shinyApp(ui = ui, server = server)









