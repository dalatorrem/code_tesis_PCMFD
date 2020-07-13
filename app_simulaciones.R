library(shiny)
library(shinythemes)


library(readxl)
library(abind)
library(rmutil)
library(dplyr)
#### paquetes para datos funcionales
library(fda)
library(funData)
library(MFPCA)


### subindices:
# l      para no. de componentes principales    d  d_no_comp_prin
# i      para tamano de la muestra              N
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
  N<-dim(scores)[1] # tamano de la muestra
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
# N<- tamano de la muestra
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

#########################################################################################################
#########################################################################################################
## Funciones para simulaciones
#########################################################################################################
#########################################################################################################
preparando_and_plot_log_val=function(N_tamano,d_divisiones,p1,p2,p3,p4,p5,p6,p7,p8,p9)
{
  t<-(0:d_divisiones)/d_divisiones
  n_len_t<-length(t)
  p1<-eval(parse(text=p1))
  p2<-eval(parse(text=p2))
  p3<-eval(parse(text=p3))
  p4<-eval(parse(text=p4))
  p5<-eval(parse(text=p5))
  p6<-eval(parse(text=p6))
  p7<-eval(parse(text=p7))
  p8<-eval(parse(text=p8))
  p9<-eval(parse(text=p9))
  mat_componentes<-matrix(c(p1,p2,p3,p4,p5,p6,p7,p8,p9),nrow=n_len_t)
  verif_arreg_de_transformaciones(mat_componentes,t)
  arreglo_de_transformaciones<-verif_arreg_de_transformaciones(mat_componentes,t)
  return(arreglo_de_transformaciones)
}
#########################################################################################################
##
## Grafico de medias
#########################################################################################################
medias_plot=function(d_divisiones,
                     fun_media_1_parte_1,fun_media_2_parte_1,fun_media_3_parte_1,
                     fun_media_1_parte_2,fun_media_2_parte_2,fun_media_3_parte_2)
{
  t<-(0:d_divisiones)/d_divisiones
  fun_media_1_parte_1<-eval(parse(text=fun_media_1_parte_1))
  fun_media_2_parte_1<-eval(parse(text=fun_media_2_parte_1))
  fun_media_3_parte_1<-eval(parse(text=fun_media_3_parte_1))
  fun_media_1_parte_2<-eval(parse(text=fun_media_1_parte_2))
  fun_media_2_parte_2<-eval(parse(text=fun_media_2_parte_2))
  fun_media_3_parte_2<-eval(parse(text=fun_media_3_parte_2))
  Conf3x2 = matrix(c(1:3), nrow=1, byrow=TRUE)
  layout(Conf3x2)
  layout.show(3)
  y_max=max(fun_media_1_parte_1,fun_media_2_parte_1,fun_media_3_parte_1,
            fun_media_1_parte_2,fun_media_2_parte_2,fun_media_3_parte_2)
  y_min=min(fun_media_1_parte_1,fun_media_2_parte_1,fun_media_3_parte_1,
            fun_media_1_parte_2,fun_media_2_parte_2,fun_media_3_parte_2)
  plot(t,fun_media_1_parte_1,"l",col=rgb(98/255,120/255,1),main="",ylab="m1",ylim=c(y_min,y_max))
  lines(t,fun_media_1_parte_2,"l",col=rgb(0,0,1))
  plot(t,fun_media_2_parte_1,"l",col=rgb(98/255,120/255,1),main="Funciones medias de los procesos",ylab="m2",ylim=c(y_min,y_max))
  lines(t,fun_media_2_parte_2,"l",col=rgb(0,0,1))
  plot(t,fun_media_3_parte_1,"l",col=rgb(98/255,120/255,1),main="",ylab="m3",ylim=c(y_min,y_max))
  lines(t,fun_media_3_parte_2,"l",col=rgb(0,0,1))
}
#########################################################################################################
##
## Funcion simular,plots,deteccion,plots
#########################################################################################################
simulacion_and_plot_con_happ=function(semilla,N_tamano,d_divisiones,p1,p2,p3,p4,p5,p6,p7,p8,p9,punto1,
                                      fun_media_1_parte_1,fun_media_2_parte_1,fun_media_3_parte_1,
                                      fun_media_1_parte_2,fun_media_2_parte_2,fun_media_3_parte_2,
                                      tipo_simul,no_bases,d1,d2,d3,an_deriv)
{
  t<-(0:d_divisiones)/d_divisiones
  n_len_t<-length(t)
  argvalsList <- list( t, t,t)
  ################################################################################
  # transformaciones
  ################################################################################
  arreglo_de_transformaciones<- preparando_and_plot_log_val(N_tamano,d_divisiones,p1,p2,p3,p4,p5,p6,p7,p8,p9)
  ################################################################################
  # ruidos
  ################################################################################
  semilla<-eval(parse(text = semilla))
  if (tipo_simul=="happ")
  {
    set.seed(semilla)
    objeto_de_Y_sin_ruido <- simMultiFunData(N = N_tamano, argvals = argvalsList,
                                             eFunType = "Fourier", eValType = "linear", M = 10,
                                             type = "split")$simData
    max_dist_Y_sin_ruido<-max(sqrt(rbind(objeto_de_Y_sin_ruido@.Data[[1]]@X,
                                         objeto_de_Y_sin_ruido@.Data[[2]]@X,
                                         objeto_de_Y_sin_ruido@.Data[[3]]@X)^2))
    desv<-max_dist_Y_sin_ruido
    set.seed(semilla)
    arreglo_de_Y_con_ruido<-addError(objeto_de_Y_sin_ruido,sd=desv)
    max_dist_Y_con_ruido<-max(sqrt(rbind(arreglo_de_Y_con_ruido@.Data[[1]]@X,
                                         arreglo_de_Y_con_ruido@.Data[[2]]@X,
                                         arreglo_de_Y_con_ruido@.Data[[3]]@X)^2))
    arreglo_de_Y_con_ruido@.Data[[1]]@X<-arreglo_de_Y_con_ruido@.Data[[1]]@X/max_dist_Y_con_ruido
    arreglo_de_Y_con_ruido@.Data[[2]]@X<-arreglo_de_Y_con_ruido@.Data[[2]]@X/max_dist_Y_con_ruido
    arreglo_de_Y_con_ruido@.Data[[3]]@X<-arreglo_de_Y_con_ruido@.Data[[3]]@X/max_dist_Y_con_ruido
    p=3 # para la apliacion p=3
    arreglo_Y_gorro<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
    arreglo_Y_gorro[,1,]<-t(arreglo_de_Y_con_ruido@.Data[[1]]@X)
    arreglo_Y_gorro[,2,]<-t(arreglo_de_Y_con_ruido@.Data[[2]]@X)
    arreglo_Y_gorro[,3,]<-t(arreglo_de_Y_con_ruido@.Data[[3]]@X)
    max_transformacion<-sqrt(max(arreglo_de_transformaciones^2))
    arreglo_Y<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
    for (j in 1:n_len_t)
    {
      arreglo_Y[j,,]<-arreglo_de_transformaciones[,,j]%*%arreglo_Y_gorro[j,,]
    }
    ## dividiendo el arreglo entre la maxima entrada para que todas las entradas 
    ## queden con errores entre -1 y 1
    arreglo_Y<-arreglo_Y/max_transformacion
  }
  else
  {
    set.seed(semilla)
    p=3 # para la apliacion p=3
    arreglo_Y_gorro<-array(rnorm(n=(n_len_t*p*N_tamano),mean=0,sd=1),c(n_len_t,p,N_tamano))
    arreglo_Y<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
    for (j in 1:n_len_t)
    {
      arreglo_Y[j,,]<-arreglo_de_transformaciones[,,j]%*%arreglo_Y_gorro[j,,]
    }
  }
  ################################################################################
  # simulaciones
  ################################################################################
  n1<-eval(parse(text=punto1))
  n2<-N_tamano-n1
  fun_media_1_parte_1<-eval(parse(text=fun_media_1_parte_1))
  fun_media_2_parte_1<-eval(parse(text=fun_media_2_parte_1))
  fun_media_3_parte_1<-eval(parse(text=fun_media_3_parte_1))
  fun_media_1_parte_2<-eval(parse(text=fun_media_1_parte_2))
  fun_media_2_parte_2<-eval(parse(text=fun_media_2_parte_2))
  fun_media_3_parte_2<-eval(parse(text=fun_media_3_parte_2))
  
  arreglo_de_fun_media_1_parte_1<-kronecker(t(rep(1,n1)),fun_media_1_parte_1)
  arreglo_de_fun_media_2_parte_1<-kronecker(t(rep(1,n1)),fun_media_2_parte_1)
  arreglo_de_fun_media_3_parte_1<-kronecker(t(rep(1,n1)),fun_media_3_parte_1)
  ##
  arreglo_de_fun_media_1_parte_2<-kronecker(t(rep(1,n2)),fun_media_1_parte_2)
  arreglo_de_fun_media_2_parte_2<-kronecker(t(rep(1,n2)),fun_media_2_parte_2)
  arreglo_de_fun_media_3_parte_2<-kronecker(t(rep(1,n2)),fun_media_3_parte_2)
  ##
  ##
  arreglo_de_fun_media_1<-cbind(arreglo_de_fun_media_1_parte_1,arreglo_de_fun_media_1_parte_2)
  arreglo_de_fun_media_2<-cbind(arreglo_de_fun_media_2_parte_1,arreglo_de_fun_media_2_parte_2)
  arreglo_de_fun_media_3<-cbind(arreglo_de_fun_media_3_parte_1,arreglo_de_fun_media_3_parte_2)
  
  datos<-array(rep(0,n_len_t*p*N_tamano),c(n_len_t,p,N_tamano))
  datos[,1,]<-arreglo_de_fun_media_1 + arreglo_Y[,1,]
  datos[,2,]<-arreglo_de_fun_media_2 + arreglo_Y[,2,]
  datos[,3,]<-arreglo_de_fun_media_3 + arreglo_Y[,3,]
  ################################################################################
  # plot ruidos
  ################################################################################
  Conf3x2 = matrix(c(1,2,3,4,5,6,7,7,7,8,9,10,11,12,12), nrow=5, byrow=TRUE)
  graf<-layout(Conf3x2)
  layout.show(graf)
  y_min=min(arreglo_Y)
  y_max=max(arreglo_Y)
  plot(t,arreglo_Y[,1,1],"l",col=gray(1),main="",ylab="Y1",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    gris=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,arreglo_Y[,1,i],"l",col=gray(gris))
  }
  plot(t,arreglo_Y[,2,1],"l",col=gray(1),main="Ruidos con dependencia",ylab="Y2",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    gris=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,arreglo_Y[,2,i],"l",col=gray(gris))
  }
  plot(t,arreglo_Y[,3,1],"l",col=gray(1),main="",ylab="Y3",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    gris=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,arreglo_Y[,3,i],"l",col=gray(gris))
  }
  ################################################################################
  # plot simulaciones
  ################################################################################
  y_max=max(datos)
  y_min=min(datos)
  plot(t,datos[,1,1],"l",col=rgb(1,0.2,0.2),main="",ylab="X1",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,1,i],"l",col=rgb(escala_color,escala_color,1))
  }
  plot(t,datos[,2,1],"l",col=rgb(1,0.2,0.2),main="Simulaciones procesos 3-variados
       (FIN DE SIMULACIoN)",ylab="X2",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,2,i],"l",col=rgb(escala_color,escala_color,1))
  }
  plot(t,datos[,3,1],"l",col=rgb(1,0.2,0.2),main="",ylab="X3",ylim=c(y_min,y_max))
  for (i in 2:N_tamano)
  {
    escala_color=0.8*(N_tamano-i)/N_tamano+0.2
    lines(t,datos[,3,i],"l",col=rgb(escala_color,escala_color,1))
  }
  ################################################################################
  # termina simulacion empieza deteccion
  ################################################################################
  lista_com_principales<-componentes_principales_multivariadas_continuas(datos)
  # valores propios
  valores_propios<-lista_com_principales$array_valores_pro
  ## plots de los valores propios estimados
  ymax<-max(log(valores_propios))
  ymin<-min(log(valores_propios))
  plot(t,log(valores_propios[1,]),"p",col=1,
       main=" INICIO DEL ANaLISIS: DETECCIoN DE PUNTO DE CAMBIO.
       Logaritmo natural  de los valores propios de B(t)B(t)' (estimados y suavizamiento) 
       (comparar con primer plot si selecciono ruidos blancos transformados)",ylim=c(ymin,ymax),ylab="")#Valores propios en el tiempo
  lines(t,log(valores_propios[2,]),"p",col=2)
  lines(t,log(valores_propios[3,]),"p",col=3)
  ## suavizamiento de los valores propios por b-splines
  ajuste_val_1<-lm(log(valores_propios[1,])~bs(t,round(d_divisiones/2)))
  ajuste_val_2<-lm(log(valores_propios[2,])~bs(t,round(d_divisiones/2)))
  ajuste_val_3<-lm(log(valores_propios[3,])~bs(t,round(d_divisiones/2)))
  
  lines(t,ajuste_val_1$fitted.values,col=1)
  lines(t,ajuste_val_2$fitted.values,col=2)
  lines(t,ajuste_val_3$fitted.values,col=3)
  # transformando los procesos
  lista_transformados<-cambio_de_base_comp_principales(datos)
  proceso_transformado<-lista_transformados$nuevo_proceso
  no_bases<-eval(parse(text = no_bases))
  # creando bases y suavizando
  Bfourier=create.fourier.basis(rangeval = c(0,1),nbasis = no_bases)
  smooth_data_1<-Data2fd(argvals = t,y=proceso_transformado[,1,],basisobj=Bfourier,dfscale = "gcv")
  smooth_data_2<-Data2fd(argvals = t,y=proceso_transformado[,2,],basisobj=Bfourier,dfscale = "gcv")
  smooth_data_3<-Data2fd(argvals = t,y=proceso_transformado[,3,],basisobj=Bfourier,dfscale = "gcv")
  # si se realiza analisis sobre la derivada.
  if (an_deriv)
  {
    smooth_data_1$coefs<-matrix_derivada(no_bases)%*%as.matrix(smooth_data_1$coefs)
    smooth_data_2$coefs<-matrix_derivada(no_bases)%*%as.matrix(smooth_data_2$coefs)
    smooth_data_3$coefs<-matrix_derivada(no_bases)%*%as.matrix(smooth_data_3$coefs)
  }
  # componentes principales funcionales
  d1<-eval(parse(text = d1))
  d2<-eval(parse(text = d2))
  d3<-eval(parse(text = d3))
  comp_principales_1<-pca.fd(smooth_data_1,nharm = d1,
                             harmfdPar = fdPar(Bfourier),centerfns = TRUE)
  comp_principales_2<-pca.fd(smooth_data_2,nharm = d2,
                             harmfdPar = fdPar(Bfourier),centerfns = TRUE)
  comp_principales_3<-pca.fd(smooth_data_3,nharm = d3,
                             harmfdPar = fdPar(Bfourier),centerfns = TRUE)
  barplot(cumsum(100*comp_principales_1$varprop), main="",xlab="No. de componentes",
          ylab="% var explicada",names.arg=1:d1,ylim=c(0,101))
  abline( h =95,col="red")
  abline( h =99,col="red")
  barplot(cumsum(100*comp_principales_2$varprop), 
          main="Estimado de la Variabilidad Explicada 
          (recomendada mayor a 95% menor a 99%)",
          xlab="No. de componentes",ylab="% var explicada",names.arg=1:d2,ylim=c(0,101))
  abline( h =95,col="red")
  abline( h =99,col="red")
  barplot(cumsum(100*comp_principales_3$varprop), main="",xlab="No. de componentes",
          ylab="% var explicada",names.arg=1:d3,ylim=c(0,101))
  abline( h =95,col="red")
  abline( h =99,col="red")
  T_estadistico_1<-estadistico_T_N_x(comp_principales_1$scores,comp_principales_1$values,FALSE)
  T_estadistico_2<-estadistico_T_N_x(comp_principales_2$scores,comp_principales_2$values,FALSE)
  T_estadistico_3<-estadistico_T_N_x(comp_principales_3$scores,comp_principales_3$values,FALSE)
  T_estadistico<-T_estadistico_1+T_estadistico_2+T_estadistico_3
  plot(seq(1,N_tamano,1),T_estadistico,
       main="Estadistico T
       (Posible punto de cambio)",
       ylab="","l",xlab="",ylim=c(0,1.7*max(T_estadistico)))
  points(which.max(T_estadistico),max(T_estadistico),pch=23,col="green",bg="red",cex = 3.5)
  coor=paste0("Candidato a PC en:",toString(which.max(T_estadistico)))
  text(x=N_tamano/2,y=1.3*max(T_estadistico),
       label= coor,cex = 1.4)
  estadistico_prueba=sum(T_estadistico)/length(T_estadistico)
  y_max=1.2*max(tabla_distri_kd,estadistico_prueba)
  d=d1+d2+d3
  x_max=max(tabla_distri_kd[,d],estadistico_prueba)
  plot(x=estadistico_prueba,y=d,pch=17,col="green",cex = 4,
       main="Distribucion Kd y Estadistico de prueba",
       ylab="No. de componentes: d1+d2+d3",xlab="valor",
       ylim=c(d-1,d+4),xlim=c(0,x_max))
  abline(h =d,col=gray(0.7))
  points(x=tabla_distri_kd[2,d],y=d,pch=4,cex=4,col=gray(0.1))
  points(x=tabla_distri_kd[3,d],y=d,pch=8,cex=4,col=gray(0.1))
  points(x=tabla_distri_kd[1,d],y=d,pch=3,cex=4,col=gray(0.1),"p")
  
  legend( 0,d+4, c("Est. Prueba",NA,"Percentil 90", "Percentil 95","Percentil 99"),
          pch = c(17 ,NA,3, 4,8),col=c("green",NA,col=gray(0.1),col=gray(0.1),col=gray(0.1)))
}
#########################################################################################################
#########################################################################################################







# Define UI for app that draws a histogram ----
ui <- fluidPage(
#######################################################################
# titulo
#######################################################################
  titlePanel("Simulaciones: Deteccion de Punto
             de Cambio en la Media para Datos Funcionales Multivariados"),
  sidebarLayout(
    sidebarPanel(style = "overflow-y:scroll;max-height: 90vh;  position:relative;",
      #######################################################################
      # Procedimientos
      #######################################################################
      helpText(h6("Se simulan N procesos IID 3-variados, (las componentes tienen dependencia),   en el intervalo [0,1]. 
      Se puede seleccionar d, el numero de divisiones  uniformes para el intervalo.
      El punto de cambio para  la media y las respectivas funciones medias. 
      La metodologia propuesta en la tesis detecta si hay punto de cambio y si lo hay,  lo estima.")),
      #######################################################################
      # Parametros
      #######################################################################
      # Input: tamano muestra ----
      h6(strong('N: Tamano de la muestra.')),
      sliderInput(inputId = "N_tamano",label = "",min = 20,max = 200,value = 50),
      # Input: numero de divisiones ----
      h6(strong('d: Numero de divisiones del intervalo.')),
      helpText(h6("Si por ejemplo selecciona d=100, se simularan los p-procesos en 0,0.01,0.02,...,0.99,1")),
      sliderInput(inputId = "d_divisiones",label = "",min = 20,max = 100,value = 100),
      h6(strong('Punto de cambio:')),
      helpText(h6("Introduzca el  punto en el que desea que cambie la media. Por ejemplo, si N=40, puede seleccionar punto=15
               si no hay punto de cambio coloque punto=N")),
      textInput(inputId = "punto1",label = "",value = 15),
      h6(strong('Tipo de simulacion')),
      helpText(h6("Si se selecciona Happ y escalada, se usaran procesos simulados como los estudiados por
         Happ y se haran escalamientos para que los ruidos queden entre 0 y 1, en este caso no se
         puede recuperar los valores propios de la matriz B(t). Si se selecciona, ruidos blancos 
         transformados, se puede recuperar aproximadamente los valores propios de la matriz B(t)")),
      selectInput(inputId = "tipo_simul",label = "",
                  choices = c("Happ y escalada"="happ","ruidos blancos transformados"="rbt")),
      #######################################################################
      # Funciones
      #######################################################################
      h6(strong("Defina la funcion B(t):")),
      helpText(h6("1ra columna, respectivamente 1ra,2da,3ra fila")),
      textInput(inputId="p1", label="", value = "t*2"),
      textInput(inputId="p2", label="", value = "t+1"),
      textInput(inputId="p3", label="", value = "t^2"),
      helpText(h6("2da columna, respectivamente 1ra,2da,3ra fila")),
      textInput(inputId="p4", label="", value = "exp(t)"),
      textInput(inputId="p5", label="", value = "sin(3*t)+2"),
      textInput(inputId="p6", label="", value = "cos(4+t)+2"),
      helpText(h6("3ra columna, respectivamente 1ra,2da,3ra fila")),
      textInput(inputId="p7", label="", value = "t^3"),
      textInput(inputId="p8", label="", value = "0*t+1.5 "),
      textInput(inputId="p9", label="", value = "t^2+1"),
      h6(strong("Media Antes del punto de cambio (Azul claro):")),
      helpText(h6("Respectivamente para 1ra,2da,3ra componente.")),
      textInput(inputId="fun_media_1_parte_1", label="", value = "t+exp(2*t)"),
      textInput(inputId="fun_media_2_parte_1", label="", value = "2*cos(3*(2*pi)*t)+sin(9*(2*pi)*t)"),
      textInput(inputId="fun_media_3_parte_1", label="", value = "2*sin((2*pi)*t)+sin(5*(2*pi)*t)+-1"),
      h6(strong("Media Despues del punto cambio:(Azul medio)")),
      helpText(h6("Respectivamente para 1ra,2da,3ra componente.")),  
      textInput(inputId="fun_media_1_parte_2", label="", value = "t+exp(2*t)+0.025*sin(12*pi*t)"),
      textInput(inputId="fun_media_2_parte_2", label="", value = "(1+0.025)*(2*cos(3*(2*pi)*t)+sin(9*(2*pi)*t))"),
      textInput(inputId="fun_media_3_parte_2", label="", value = "2*sin((2*pi)*t)+sin(5*(2*pi)*t)-1+0.025"),
      #######################################################################
      # Semilla
      #######################################################################
      h6(strong("Semilla para simulaciones:")),
      textInput(inputId="semilla", label="", value = 1234),
      h5(strong("DETECCIoN DE PUNTO DE CAMBIO")),
      h6(strong("Numero de Funciones de la base:")),
      helpText("En esta apliaccion se esta trabajando con bases de Fourier. 
               El parametro de regularizacion es obtenido por GCV."),
      textInput(inputId="no_bases", label="", value = 15),
      helpText(h6("Esta aplicacion permite realizar un analisis usando las derivadas teniendo en cuenta que se trabaja en la base de Fourier.
               Este analisis permite detectar puntos de cambio en la forma de la media.")),
      checkboxInput(inputId="an_deriv", label="Analizar derivadas", value = FALSE),
      h6(strong("Numero de componetes:")),
      helpText(h6("Indique el numero de componentes principales para el ACP univariado funcional.
                Debe ser menor que el numero de bases. Escoja un numero de componentes que logre
                explicar entre el 95% y el 99%. En esta aplicacion se puede decidir si la suma de 
                los tres valores es menor o igual a 30. Respectivamente para 1ra,2da,3ra componente. 
               Tenga en cuenta el grafico: Estimado de la Variabilidad Explicada.")),
      textInput(inputId="d1", label="", value =10),
      textInput(inputId="d2", label="", value =10),
      textInput(inputId="d3", label="", value =10)
    ),
    mainPanel(style = "overflow-y:scroll;max-height: 90vh;  position:relative;",
      fluidRow(
        plotOutput(outputId="plot_log_val_prop",height = "220px"),  
        plotOutput(outputId="plot_medias",height = "200px"),
        plotOutput(outputId="plot_simulaciones",height = "900px")
      )
    )
      ))

# Define server ----
server <- function(input, output) {
  output$plot_log_val_prop <- renderPlot({
    preparando_and_plot_log_val(input$N_tamano,input$d_divisiones,input$p1,input$p2,input$p3,input$p4,input$p5,input$p6,input$p7,
                                input$p8,input$p9)
  })
  output$plot_medias <- renderPlot({
    medias_plot(input$d_divisiones,
                input$fun_media_1_parte_1,input$fun_media_2_parte_1,input$fun_media_3_parte_1,
                input$fun_media_1_parte_2,input$fun_media_2_parte_2,input$fun_media_3_parte_2)
    })
  output$plot_simulaciones <- renderPlot({
    simulacion_and_plot_con_happ(input$semilla,input$N_tamano,input$d_divisiones,input$p1,input$p2,input$p3,
                        input$p4,input$p5,input$p6,input$p7,input$p8,input$p9,input$punto1,
                        input$fun_media_1_parte_1,input$fun_media_2_parte_1,input$fun_media_3_parte_1,
                        input$fun_media_1_parte_2,input$fun_media_2_parte_2,input$fun_media_3_parte_2,
                        input$tipo_simul,input$no_bases,input$d1,input$d2,input$d3,input$an_deriv)
  })
}
shinyApp(ui = ui, server = server)
