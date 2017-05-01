#Programa para estudiar un conjunto de osciladores harmónicos, acoplados por un potencial por pares
#En 1D
#Este programa calcula las interacciones tipo Lennard Jones con todos los osciladores de la red
#Tiene la posibilidad de modificar los parámetros de la interacción para todos los elementos de la red, de forma separada se puede cambiar el parámetro de la parte repulsiva y de la parte atractiva
#La masa de cada oscilador se pueden especificar
#Las constantes de oscilador se pueden también especificar de forma separada. La forma en la cual se cuentan los osciladores es:
# m_1 ---- k_1 ----- m_2 ----- k_2 ---- m_3 ... m_(N-2) ---- k_(N-2) ---- m_(N-1) ---- k_(N-1) ---- m_N

#Encabezado no es necesario modificar
#Librerías
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
using DifferentialEquations
using Plots
using Distributions
#Generador de semilla de los números aleatorios
srand(123)
#/Encabezado

#Sistema
N = 4; #Tamaño del Sistema
a0 = 10.; #Parámetro de la red
tiempo = (0,10.); #Intervalo donde se calculará la solución
k = Array{Float64}(N-1); #Vector de parámetros del potencial de oscilador. Para llenarlo con una distribución alrededor de un número, usa la función D
# k = 1.4 # Para el caso de sólo un parámetro. Se suponen condiciones a la frontera libres
#a = Array{Float64}(N,N); #Parámetros de la parte no lineal del potencial, caso interacciones diferentes
#b = Array{Float64}(N,N); #Parámetros opcionales para otra cosa, casi interacciones diferentes
a = 0.001; #Amplitud de la interacción Lennard-Jones
b = 0.001; #Mínimo de la interacción Lennard-Jones
x = Array{Float64}(N); #Posiciones canónicas del oscilador
v = Array{Float64}(N); #Velocidades canónicas del oscilador
xo = Array{Float64}(N); #Posiciones iniciales
vo = Array{Float64}(N); #Velocidades iniciales
nl = zeros{Float64}(2*N); #Interacción no lineal
m = Array{Float64}(N); #Masas
#/Sistema

#Opciones para las fluctuaciones en k,a y b:
kµ = 10. #Promedio k
kσ = 1.  #Desvest k o parámetro de Lorenziana k
#aµ = 10. #Promedio a
#aσ = 1.  #Desvest a o parámetro de Lorenziana a
#bµ = 10. #Promedio b
#bσ = 1.  #Desvest b o parámetro de Lorenziana b
#/Opciones para las fluctuaciones en k,a y b

#Funciones con parámetros que puedes cambiar
#Pontencial
function LJ(x1,x2,a,b,a0) #Fuerza Lennard-Jones, x1 y x2 se sugieren usar como variables y a y b como parámetros de la interacción. a0 es el parámetro de la red
  r=(a0-x1+x2)
  return 12.*a*((b/r)^12-(b/r)^6)/r #a mide la magnitud de la fuerza mientras que b localiza espacialmente el mínimo del potencial
end
#/P
#Distribuciones
function D(N,µ,σ,β) #Distribución puede ser usada para generar k's que varíen un poco entre ellas
  dist = Normal(µ,σ) #Gaussiana
#  dist = Cauchy(µ,β) #Lorentz
  k = rand(dist,N)
  return k
end
#/D
#/Funciones con parámetros que puedes cambiar

#Definiciones necesarias para el programa. Matriz A, condiciones iniciales,
k = D(N-1,kµ,kσ,0)

#Generación de la matriz del sistema. Tamaño: 2N*2N, incluye el par de ecuaciones de Hamilton para cada oscilador
X = zeros{Float64}(2*N) #El vector X contiene las posiciones en las primeras N entradas y las velocidades en el resto
A = zeros{Float64}(2*N,2*N) #La matriz A permite resolver simultáneamente las ecuaciones de Hamilton
#/G
#Poniendo las condiciones iniciales en X
for i in (1:N)
  X[i] = xo[i]
  X[N+i] = vo[i]
end
#/P
#Construcción de la matriz A y definición del problema. Ver A.pdf

#Ecuaciones Dx = v
for i in (1:N)
  A[i,N+i]= 1
end
#/Ecuaciones Dx = v

#Ecuaciones Dv = F/m
#Frontera
A[N+1,1] = -k[1]/m[1]
A[N+1,2] = k[1]/m[1]
A[2*N,2*N-1] = k[N-1]/m[N]
A[2*N,2*N] = -k[N-1]/m[N]
#/Frontera
#Interior de la cadena
for i in (N+2:2*N-1)
  A[i,i-1] = k[i-1]/m[i-N]
  A[i,i] = -(k[i]+k[i-1])/m[i-N]
  A[i,i+1] = k[i]/m[i-N]
end
#/Interior

#ECUACIÓN DIFERENCIAL:
F(t,u) = A*u + f(u) #f(u) es la parte no lineal EN EL FUTURO, PROBAR INTRODUCIR a0 como parámetro para ver si es más rápido
#parte no lineal de la ecuación diferencial
function f(u) #sólo vamos a usar las posiciones
  for i in (1:N)
    nl[i+N]=0.
    for j in (i+1:N) #Es lento en el sentido de que calcula dos veces U(x,y) y U(y,x), que cuando las interacciones son simétricas pues son iguales PERO no usa memoria adicional
      #nl[i+N]+=U(i,j,a[i,j],b[i,j],a0) #Interacciones diferentes
      nl[i+N]+=U(i,j,a,b,a0) #Interacciones iguales
    end
    for j in (1:i-1)
      #nl[i+N]+=U(i,j,a[i,j],b[i,j],a0) #Interacciones diferentes
      nl[i+N]+=U(i,j,a,b,a0)
    end
  end
  return nl
end
#/ECUACIÓN DIFERENCIAL

#PROGRAMA DE JULIA
problema = ODEProblem(F,X,tiempo)
solucion = solve(problema)

#Análisis de la solución
E = Array{Float64}(N) #Energía en los modos
ω = Array{Float64}(N) #Frecuencias normales
base = Array{Float64}(N) #Modos normales

for i in (1:N)
  ω[i] = 2*sin(i*pi/2(N+1))
end

for i in (1:N)
  Ax = 0;
  Av = 0;
  for j in ()
  E[i] = 0.5*
