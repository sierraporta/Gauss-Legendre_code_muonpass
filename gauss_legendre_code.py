# -*- coding: utf-8 -*-
# /************************************************************************/
# /*                                                                      */
# /* Script:  gauss_legendre_code.py                                      */
# /* Issue:   Resolviendo una integral por el método de Gauss-Legendre    */
# /*                                                                      */
# /************************************************************************/
# /* Authors:  David Sierra Porta                                         */
# /* e-mail:   dsierrap@uis.edu.co                                        */
# /*                                                                      */
# /************************************************************************/
# /* Comments: Resolviendo integral numéricamente                         */
# /*                                                                      */
# /************************************************************************/

from numpy import *
 
##################################################################
# Vamos a definir de manera recursiva todos los polinomios de orden n
def Legendre(n,x):
 x=array(x)
 if (n==0):
  return x*0+1.0
 elif (n==1):
  return x
 else:
  return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n
 
##################################################################
# Definimos las derivadas d elos polinomios de Legendre en terminos de los mismos polinomios
def DLegendre(n,x):
 x=array(x)
 if (n==0):
  return x*0
 elif (n==1):
  return x*0+1.0
 else:
  return (n/(x**2-1.0))*(x*Legendre(n,x)-Legendre(n-1,x))
##################################################################
# Raices de los polinnomios obtenidas usando metodo de Newton-Raphson
def LegendreRoots(polyorder,tolerance=1e-20):
 if polyorder<2:
  err=1 # No se puede obtener raiz
 else:
  roots=[]
# Los polinomios son funciones pares y impares alternadas. Evaluamos en un numero semientero de raices. 
  for i in range(1,int(polyorder)/2 +1):
   x=cos(pi*(i-0.25)/(polyorder+0.5))
   error=10*tolerance
   iters=0
   while (error>tolerance) and (iters<1000):
    dx=-Legendre(polyorder,x)/DLegendre(polyorder,x)
    x=x+dx
    iters=iters+1
    error=abs(dx)
   roots.append(x)
# Uso de simetrias para obtener otras raices
  roots=array(roots)
  if polyorder%2==0:
   roots=concatenate( (-1.0*roots, roots[::-1]) )
  else:
   roots=concatenate( (-1.0*roots, [0.0], roots[::-1]) )
  err=0 # Raices determinadas correctamente
 return [roots, err]
##################################################################
# Coeficientes de los pesos para la expansion
def GaussLegendreWeights(polyorder):
 W=[]
 [xis,err]=LegendreRoots(polyorder)
 if err==0:
  W=2.0/( (1.0-xis**2)*(DLegendre(polyorder,xis)**2) )
  err=0
 else:
  err=1 # El error es muy grande no pueden determinarse las raices, tampoco los pesos
 return [W, xis, err]
##################################################################
# La integral a resolver 
# func 	: el integrando
# a, b 	: limites inferior y superior de la integral
# polyorder 	: orden el poinomio de Legendre a usar
#
def GaussLegendreQuadrature(func,polyorder,a,b):
 [Ws,xs, err]= GaussLegendreWeights(polyorder)
 if err==0:
  ans=(b-a)*0.5*sum( Ws*func( (b-a)*0.5*xs+ (b+a)*0.5 ) )
 else: 
  # (en caso de error)
  err=1
  ans=None
 return [ans,err]
##################################################################
# Integrando, lo cambiamos a conveniencia
# Se definen tambien aca los parametros de la funcion a integrar
# Aca ponemos los limites de la integral explicitamente (a,b)
alpha=2.3E+00	
beta=2.57E-06
a=3
b=6
def func(x):
 return 1/(alpha+beta*x)
##################################################################
# Definimos el orden de aproximacion
order=5
[Ws,xs,err]=GaussLegendreWeights(order)
if err==0:
 print "Orden    : ", order
 print "Raices   : ", xs
 print "Pesos    : ", Ws
else:
 print "No se pueden evaluar las raices/pesos"
 
# Integrando la funcion
[ans,err]=GaussLegendreQuadrature(func,order,a,b)
if err==0:
 print "Integral : ", ans
else:
 print "No se puede evaluar la integral"
