library(readxl)
library(abind)
library(rmutil)
library(dplyr)
#### paquetes para datos funcionales
library(fda)
#library(funData)
library(MFPCA)


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
## Estadístico S_d ver Berkes et al Detecting changes in the mean of functional observations pag 932, 3.4
#########################################################################################################
#########################################################################################################

# ingresa la matriz de scores después de haber realizado componentes principales para datos 
# funcionales, el obejto que entra tiene la estructura de PCA$scores, se debe tener cuidado 
# con la estructura si se usa otra función de componentes principales funcionales

estadistico_S_N_D=function(scores,valores_propios)
{
  N<-dim(scores)[1] # tamaño de la muestra
  d_no_comp_prin<-dim(scores)[2]
  # prueba función 1
  s=0
  for (l in 1:d_no_comp_prin)
  {
    inv_val<-1/valores_propios[l]
    score_l<-scores[,l]
    suma=0
    for (i in 1:N)
    {
      sub_1<-score_l[1:i]
      a<-sum(sub_1) # primera suma dentro del parentesis que está al cuadrado
      b<-sum(score_l) # segunda suma dentro del parentesis que está al cuadrado
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
## Estadístico T_d ver Berkes et al Detecting changes in the mean of functional observations pag 931, 3.3
#########################################################################################################
#########################################################################################################

# ingresa la matriz de scores después de haber realizado componentes principales para datos 
# funcionales, el obejto que entra tiene la estructura de PCA$scores, se debe tener cuidado 
# con la estructura si se usa otra función de componentes principales funcionales

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
      a<-sum(sub_1) # primera suma dentro del parentesis que está al cuadrado
      b<-sum(score_l) # segunda suma dentro del parentesis que está al cuadrado
      parentesis_al_cuadrado<-(a-(i/N)*b)^2
      suma=suma+inv_val*parentesis_al_cuadrado
    }
    T_i<-suma/N
    est_t_x<-rbind(est_t_x,T_i)
  }
  if (plot_bool==TRUE)
  {
  plot(seq(1,N,1),est_t_x,main="Gráfico estadístico T","l")
  }
  return(est_t_x)
}   



#########################################################################################################
#########################################################################################################
##  SE INGRESA UN VALOR K ENTERO Y DEVUELVE LA MATRIZ DE DERIVADA PARA UNA BASE DE FOURIER DE 
##  ORDEN K, LA BASE DE FOURIER ESTÁ EN EL ORDEN DADO POR LA FUNCIÓN create.fourier.basis()$names
#########################################################################################################
#########################################################################################################



matrix_derivada<-function(no_bases)
{
  ceros<-t(as.matrix(rep(0,no_bases)))
  mat_derivada<-ceros
  for (n_b in 2:no_bases) 
  {
    # para las posiciones pares o sea para los senos lo derivada es la contracción del seno (2pi)
    # y el signo es positivo, para el coseno es análogo salvo que el signo es negativo. para 1 der=0
    der<-floor(n_b/2)*2*pi*(-1)^n_b  
    derivada<-ceros
    # para las posiciones pares (senos) lo envía uno adelante o sea al coseno del mismo argumento
    # para las posiciones impares (cosenos) lo envía uno atras o sea al seno del mismo argumento
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
# las medidas tienen que ajustarse al intervalo [0,1] haciendo una traslación y un reescalamiento
# del tiempo. 
# en esta función se calculan los componentes principales p-variados en cada división de tiempo
# del interavlo [0,1]
# el argumento procesos es un arreglo [len,p,N]
# len<- número de intervalos en [0,1]
# p<- dimensiones del proceso (p-variado)
# N<- tamaño de la muestra
# la función retorna una lista en la que hay:
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
      # los vectores propios lo más cercano a esa base, dado que hay 2^p posibilidades para seleccionar 
      # la base ortogonal de vectores propios, en el paso 1 se selecciona como los vectores canónicos
    }else
    {
      vectores_propios_t_ant<-vectores_propios_t # para evitar al máximo cambios bruscos en 
      # la matriz ortonormal con la base de vectores propios se selecciona la base de vectores propios 
      # que esté lo más cerca posible a la anterior base 
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
#### Seleeción de la base de vectores propios para garantizar continuidad
#########################################################################################################
#########################################################################################################

# esta función recibe dos matrices que serian la base ortonormal en el instante de tiempo anterior y la 
# de t. La función hace todas las posibles de los signos de las columnas y devuelve aquella matrix que 
# tenga la menor distancia de frobenius a la base ortonormal en el instante anterior.


seleccion_vectores_propios_continuidad<-function(vectores_propios_t_ant,vectores_propios_t_des)
{
  p<-dim(vectores_propios_t_des)[1]  # dimensión
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

# Esta función recibe los procesos como un arreglo de la misma forma que las demás funciones
# En la función se realiza en cada instante de tiempo una transformación que deje la matriz 
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
#### función para verificar si las componentes de una matriz cumplen las condiciones de Berrendero  
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
  # deben cumplir las hipótesis de Berrendero, o sea los lambda en el tiempo deben ser continuos, 
  # deben ser una sucesión decreciente de funciones
  
  valores_propios<-matrix(rep(0,p*n_len_t),ncol = n_len_t)
  for (j in 1:n_len_t)
  {
    valores_propios[,j]<-eigen(arreglo_de_covarianzas[,,j])$values
  }
  
  
  # gráficos de las  funciones de los valores propios
  
  
  #ymax<-max((valores_propios))
  #ymin<-min((valores_propios))
  
  #plot(t,(valores_propios[1,]),ylim = c(ymin,ymax),"l",col=1, 
  #     main="",ylab="")# gráfico de valores propios
  #for (k in 2:p)
  #{
  #  lines(t,valores_propios[k,],col=k)
  #}

  
  # gráficos de los logaritmos de las  funciones de los valores propios
  
  
  ymax<-max(log(valores_propios))
  ymin<-min(log(valores_propios))
  
  plot(t,log(valores_propios[1,]),ylim = c(ymin,ymax),"l",col=1,
       main="Logaritmo natural  de los valores propios de B(t)B(t)' ",ylab="") # Gráfico de logaritmos naturales de los valores propios en el intervalo
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


