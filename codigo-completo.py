# -*- coding: utf-8 -*-
""" 
MAP5725 - Prof. Alexandre Roma
Aluna: Karina Yukimi Peixoto Sakurai - 9796634
Tarefa 2
"""
import math
import matplotlib.pyplot as plt
import numpy as np

def phi_Euler(t,y,dt,f):                            # Função de discretização do Método de Euler
    
    k1 = f(t, y)
    
    return k1

def phi_EulerImplicito(t,y,dt,f):                   # Função de discretização do Método de Euler Implícito
    
    A = np.array([[0,1],[-4*(t**2), 1/t]])
    I = np.array([[1,0],[0, 1]])
    M = I - dt*A
    k1 = f(t + dt, np.linalg.inv(M).dot(y))
    
    return k1

def phi_Trapezio(t,y,dt,f):                         # Função de discretização do Método do Trapézio
    
    k1 = f(t, y)
    A = np.array([[0,1],[-4*(t**2), 1/t]])
    I = np.array([[1,0],[0, 1]])
    M = I - (dt/2)*A
    N = I + (dt/2)*A
    k2 = f(t+dt, np.linalg.inv(M).dot(N).dot(y))
    
    return (k1 + k2)/2

def f(t, y):
   
    f0 = y[1]
    f1 = y[1]/t - 4*(t**2)*y[0]
    
    return np.array([f0,f1])

def sol_exata(t):                                               # Função que calcula a solução exata para um t específico

    sol = np.array([math.sin(t**2), 2*t*math.cos(t**2)])
    
    return sol

#Lendo comandos
print("Tarefa 2 \n")
print("Opções disponíveis: \n")

print("** Métodos Numéricos: **\n")

print("-> Digite '1' para Método de Euler")
print("-> Digite '2' para Método de Euler Implícito")
print("-> Digite '3' para Método do Trapézio")
print("-> Digite '4' para gerar um gráfico comparativo de métodos e solução exata \n")

print("** Tabela de convergência e erro de discretização: **\n")

print("-> Digite '5' para gerar tabela para o Método de Euler")
print("-> Digite '6' para gerar tabela para o Método de Euler Implícito")
print("-> Digite '7' para gerar tabela para o Método do Trapézio \n")

print("** Gráficos com aproximações numéricas para diferentes valores de dt: **\n")
print("-> Digite '8' para gerar aproximações usando o Método de Euler")
print("-> Digite '9' para gerar aproximações usando o Método de Euler Implícito")
print("-> Digite '10' para gerar aproximações usando o Método do Trapézio")

print("** CASO 2: Problema de Cauchy Genérico\n")

print("-> Digite '11' para gerar tabela para o Método de Euler")
print("-> Digite '12' para gerar tabela para o Método de Euler Implícito")
print("-> Digite '13' para gerar tabela para o Método do Trapézio \n")

opcao = input()

t_n = [math.sqrt(math.pi)]; T = 10;                                             # Intervalo de tempo
y_n = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Condições iniciais (x0, y0)

#Estes vetores serão utilizados caso a opção '4' for selecionada
y_n1 = [np.array([0,-2*math.sqrt(math.pi)])]
y_n2 = [np.array([0,-2*math.sqrt(math.pi)])]
y_n3 = [np.array([0,-2*math.sqrt(math.pi)])]
y_exata = [np.array([0,-2*math.sqrt(math.pi)])]

n = 2**17                                                                       # Discretização do intervalo de tempo
dt = (T-t_n[-1])/n

if (opcao != '1' and opcao != '2' and opcao != '3' and opcao != '4' and opcao != '5' and opcao != '6' and opcao != '7'
    and opcao != '8' and opcao != '9' and opcao != '10' and opcao != '11' and opcao != '12' and opcao != '13'):
    
    print("Entrada inválida. Fim do programa.")
    
else:
    
    opcao = int(opcao)

    # *********** Opções 1, 2, 3 e 4 ***********

    if (opcao <= 4): 
    
        while t_n[-1] < T:
                
            y = 0
                
            if (opcao == 1):
                    
                y = y_n[-1] + dt*phi_Euler(t_n[-1],y_n[-1],dt,f)                            #Chamada Método de Euler
                y_n.append(y)
                    
            elif (opcao == 2):
                    
                y = y_n[-1] + dt*phi_EulerImplicito(t_n[-1],y_n[-1],dt,f)                   #Chamada Método de Euler Implícito
                y_n.append(y)
                    
            elif (opcao == 3):
                    
                y = y_n[-1] + dt*phi_Trapezio(t_n[-1],y_n[-1],dt,f)                         #Chamada Método do Trapézio
                y_n.append(y)                 

            elif (opcao == 4):
                #Neste caso foi escolhido o método 4, os 3 métodos serão gerados simultaneamente. 
                y_n1.append(y_n1[-1] + dt*phi_Euler(t_n[-1],y_n1[-1],dt,f))
                y_n2.append(y_n2[-1] + dt*phi_EulerImplicito(t_n[-1],y_n2[-1],dt,f))
                y_n3.append(y_n3[-1] + dt*phi_Trapezio(t_n[-1],y_n3[-1],dt,f))
                y_exata.append(sol_exata(t_n[-1]))

            t_n.append(t_n[-1] + dt)
            dt = min(dt, T-t_n[-1])

            if (t_n[-1] == T):
                aprox = y

        if (opcao != 4): # Opções 1, 2 e 3

            y_n = np.array(y_n)
        
            # Plotando gráfico de x em função de t
            plt.plot(t_n, y_n[:,0], 'k-.', markersize=0.1)
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('x(t)')
            if (opcao == 1):
                plt.title('Método de Euler - Gráfixo x(t)')
            elif (opcao == 2):
                plt.title('Método de Euler Implícito - Gráfico x(t)')
            else:
                plt.title('Método do Trapézio - Gráfixo x(t)')
            plt.show()

            # Plotando gráfico de y em função de t
            plt.plot(t_n, y_n[:,1], 'k-.', markersize=0.1)
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('y(t)')
            if (opcao == 1):
                plt.title('Método de Euler - Gráfixo y(t)')
            elif (opcao == 2):
                plt.title('Método de Euler Implícito - Gráfico y(t)')
            else:
                plt.title('Método do Trapézio - Gráfixo y(t)')
            plt.show()

        else: # Opção 4

            y_n1 = np.array(y_n1)
            y_n2 = np.array(y_n2)
            y_n3 = np.array(y_n3)
            y_exata = np.array(y_exata)

            # Plotando gráfico de x em função de t
            plt.plot(t_n, y_n1[:,0], 'k--', markersize=0.1, label = 'Euler')
            plt.plot(t_n, y_n2[:,0], 'k-.', markersize=0.1, label = 'Euler Implícito')
            plt.plot(t_n, y_n3[:,0], 'k+', markersize=0.1, label = 'Trapézio')
            plt.plot(t_n, y_exata[:,0], 'k-', markersize=0.1, label = 'Solução Exata')
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('x(t)')
            plt.title('Comparativos de Métodos n = 2**17 - Gráfixo x(t)')
            plt.show()

            # Plotando gráfico de y em função de t
            plt.plot(t_n, y_n1[:,1], 'k--', markersize=0.1, label = 'Euler')
            plt.plot(t_n, y_n2[:,1], 'k-.', markersize=0.1, label = 'Euler Implícito')
            plt.plot(t_n, y_n3[:,1], 'k+', markersize=0.1, label = 'Trapézio')
            plt.plot(t_n, y_exata[:,1], 'k-', markersize=0.1, label = 'Solução Exata')
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('y(t)')
            plt.title('Comparativos de Métodos n = 2**17 - Gráfixo y(t)')
            plt.show()

    # *********** Opções 5, 6 e 7: Gerar tabelas de convergência e erro de discretização global ***********    

    if (opcao >= 5 and opcao <= 7):

        potencia = [2,5,10,15,18]
        linha = []
        erro = 0
        erro2 = 0
        
        for item in potencia:
            n = 2**item
            t_n = [math.sqrt(math.pi)]; T = 10;                                             # Intervalo de tempo
            y_n = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Condições iniciais (x0, y0)
            y_exata = [np.array([0,-2*math.sqrt(math.pi)])]
            dt = (T-t_n[-1])/n
            
            while t_n[-1] < T:
                if (opcao == 5):
                        
                    y_n.append(y_n[-1] + dt*phi_Euler(t_n[-1],y_n[-1],dt,f))                        # Cálculo de y usando Euler
                        
                elif (opcao == 6):
                    
                    y_n.append(y_n[-1] + dt*phi_EulerImplicito(t_n[-1],y_n[-1],dt,f))               # Cálculo de y usando Euler Implícito
                        
                elif (opcao == 7):
                        
                    y_n.append(y_n[-1] + dt*phi_Trapezio(t_n[-1],y_n[-1],dt,f))                     # Cálculo de y usando o método do trapezio                    

                y_exata.append(sol_exata(t_n[-1]))
                t_n.append(t_n[-1] + dt)
                dt = min(dt, T-t_n[-1])

            # Para cálculo da ordem de convergência
            t_n2 = [math.sqrt(math.pi)]; T = 10;                                             # Intervalo de tempo
            y_n2 = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Condições iniciais (x0, y0)
            dt = (T-t_n2[-1])/n
            dt = 2*dt

            while t_n2[-1] < T:
                if (opcao == 5):
                        
                    y_n2.append(y_n2[-1] + dt*phi_Euler(t_n2[-1],y_n2[-1],dt,f))                        # Cálculo de y usando Euler
                        
                elif (opcao == 6):
                    
                    y_n2.append(y_n2[-1] + dt*phi_EulerImplicito(t_n2[-1],y_n2[-1],dt,f))               # Cálculo de y usando Euler Implícito
                        
                elif (opcao == 7):
                        
                    y_n2.append(y_n2[-1] + dt*phi_Trapezio(t_n2[-1],y_n2[-1],dt,f))                     # Cálculo de y usando o método do trapezio
                    
                t_n2.append(t_n2[-1] + dt)
                dt = min(dt, T-t_n2[-1])
            
            erro = y_exata[len(y_exata) - 1] - y_n[len(y_n) - 1]
            erro2 = y_n2[len(y_n2) - 1] - y_n[len(y_n) - 1]
            q = erro2/erro
            linha.append([n, dt, "{:.2e}".format(abs(erro[0])), "{:.2e}".format(abs(erro2[1])), "{:.2e}".format(abs(q[0])), "{:.2e}".format(abs(q[1]))])
            
        linha = np.array(linha)
        print("n ----- dt ----- erro_x ---- erro_y ---- q_x ---- q_y")
        print(linha)
        
    # *********** Métodos 8, 9 e 10 ***********

    if (opcao >= 8 and opcao <= 10):

        potencia = [12, 15, 18]                                                             # Potências utilizadas no cálculo de n (número de passos)
        dados = []                                                                          # Variável que vai armazenar listas de y para cada n
        estilo = ["k:","k-.","k--"]                                                         # Estilo de linhas para os gráficos 

        _y = []
        _t = []
        
        for item in potencia:

            t_n = [math.sqrt(math.pi)]; T = 10;                                             # Inicializando intervalo de tempo
            y_n = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Inicializando condições iniciais (x0, y0)
            n = 2**item                                                                     # Número de passos de iteração (discretização do intervalo)
            dt = (T-t_n[-1])/n                                                              # Tamanho do passo de iteração
            y = 0                                                                           # Inicialização da variável y
                
            while t_n[-1] < T:
                    
                if (opcao == 8):
                        
                    y = y_n[-1] + dt*phi_Euler(t_n[-1],y_n[-1],dt,f)                        # Chamada método Euler para cálculo de y_{n+1}
                        
                elif (opcao == 9):
                        
                    y = y_n[-1] + dt*phi_EulerImplicito(t_n[-1],y_n[-1],dt,f)               # Chamada método Euler Implícito para cálculo de y_{n+1}
                        
                elif (opcao == 10):
                        
                    y = y_n[-1] + dt*phi_Trapezio(t_n[-1],y_n[-1],dt,f)                     # Chamada método do trapézio para cálculo de y_{n+1}

                y_n.append(y)                                                               # Adicionando y calculado na lista de valores de y
                t_n.append(t_n[-1] + dt)                                                    # Adicionando t atual na lista de valores de t
                dt = min(dt, T-t_n[-1])                                                     

            _y.append(y_n)
            _t.append(t_n)

        legenda = ""

        # Plotando gráfico de x em função de t
        for item in range(len(potencia)):

            y_n = _y[item]
            t_n = _t[item]
            y_n = np.array(y_n)
            pot = str(potencia[item])
            legenda = "n = 2**"+pot

            plt.plot(t_n, y_n[:,0], estilo[item], label = legenda, markersize=0.1) 
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('x(t)')
            plt.legend(loc='upper left')
            if (opcao == 8):
                legenda = "Euler"
                plt.title('Método de Euler - Gráfixo x(t)')
            elif (opcao == 9):
                legenda = "Euler Implícito"
                plt.title('Método de Euler Implícito - Gráfico x(t)')
            else:
                legenda = "Trapézio"
                plt.title('Método do Trapézio - Gráfixo x(t)')
            
        plt.show()

        # Plotando gráfico de x em função de t
        for item in range(len(potencia)):

            y_n = _y[item]
            t_n = _t[item]
            y_n = np.array(y_n)
            pot = str(potencia[item])
            legenda = "n = 2**"+pot

            plt.plot(t_n, y_n[:,1], estilo[item], label = legenda, markersize=0.1) 
            plt.xlabel('t (em unidades de tempo)')
            plt.ylabel('y(t)')
            plt.legend(loc='upper left')
            if (opcao == 8):
                legenda = "Euler"
                plt.title('Método de Euler - Gráfixo y(t)')
            elif (opcao == 9):
                legenda = "Euler Implícito"
                plt.title('Método de Euler Implícito - Gráfico y(t)')
            else:
                legenda = "Trapézio"
                plt.title('Método do Trapézio - Gráfixo y(t)')
            
        plt.show()
        
# *********** Opções 11, 12 e 13: CASO 2 - Problema de Cauchy Genérico ***********    

    if (opcao >= 11 and opcao <= 13):

        potencia = [2,5,10,15,18]
        linha = []
        
        for item in potencia:
            n = 2**item
            t_n = [math.sqrt(math.pi)];                                                     # Intervalo de tempo
            y_n = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Condições iniciais (x0, y0)
            y_exata = [np.array([0,-2*math.sqrt(math.pi)])]
            dt = (T-t_n[-1])/n
            
            while t_n[-1] < T:
                if (opcao == 11):
                        
                    y_n.append(y_n[-1] + dt*phi_Euler(t_n[-1],y_n[-1],dt,f))                        # Cálculo de y usando Euler
                        
                elif (opcao == 12):
                    
                    y_n.append(y_n[-1] + dt*phi_EulerImplicito(t_n[-1],y_n[-1],dt,f))               # Cálculo de y usando Euler Implícito
                        
                elif (opcao == 13):
                        
                    y_n.append(y_n[-1] + dt*phi_Trapezio(t_n[-1],y_n[-1],dt,f))                     # Cálculo de y usando o método do trapezio                    

                t_n.append(t_n[-1] + dt)
                dt = min(dt, T-t_n[-1])
            # Para cálculo da ordem de convergência
            t_n2 = [math.sqrt(math.pi)];                                                        # Intervalo de tempo
            y_n2 = [np.array([0,-2*math.sqrt(math.pi)])]                                        # Condições iniciais (x0, y0)
            t_nMeio = [math.sqrt(math.pi)];                                                     # Intervalo de tempo
            y_nMeio = [np.array([0,-2*math.sqrt(math.pi)])]                                     # Condições iniciais (x0, y0)
            dt2 = 2*(T-t_n2[-1])/n
            dtMeio = (T-t_n2[-1])/(2*n)
            while t_n2[-1] < T:
                if (opcao == 11):
                        
                    y_n2.append(y_n2[-1] + dt2*phi_Euler(t_n2[-1],y_n2[-1],dt2,f))                        # Cálculo de y usando Euler
                        
                elif (opcao == 12):
                    
                    y_n2.append(y_n2[-1] + dt2*phi_EulerImplicito(t_n2[-1],y_n2[-1],dt2,f))               # Cálculo de y usando Euler Implícito
                    
                elif (opcao == 13):
                        
                    y_n2.append(y_n2[-1] + dt2*phi_Trapezio(t_n2[-1],y_n2[-1],dt2,f))                     # Cálculo de y usando o método do trapezio
                    
                t_n2.append(t_n2[-1] + dt2)
                dt2 = min(dt2, T-t_n2[-1])

            while t_nMeio[-1] < T:
                if (opcao == 11):

                    y_nMeio.append(y_nMeio[-1] + dtMeio*phi_Euler(t_nMeio[-1],y_nMeio[-1],dtMeio,f))
                        
                elif (opcao == 12):

                    y_nMeio.append(y_nMeio[-1] + dtMeio*phi_EulerImplicito(t_nMeio[-1],y_nMeio[-1],dtMeio,f))
                    
                elif(opcao == 13):

                    y_nMeio.append(y_nMeio[-1] + dtMeio*phi_Trapezio(t_nMeio[-1],y_nMeio[-1],dtMeio,f))
                    
                t_nMeio.append(t_nMeio[-1] + dtMeio)
                dtMeio = min(dtMeio, T-t_nMeio[-1])


            numerador = y_n2[len(y_n2) - 1] - y_n[len(y_n) - 1]
            denominador = y_n[len(y_n) - 1] - y_nMeio[len(y_nMeio) - 1]
            
            resx = abs(numerador[0]/denominador[0])
            px = math.log(resx, 2)                                      #Cálculo da ordem de convergência aprox na variável x
            
            resy = abs(numerador[1]/denominador[1])
            py = math.log(resy, 2)                                      #Cálculo da ordem de convergência aprox na variavel y
            
            linha.append([item, n, dt, "{:.2e}".format(y_n[len(y_n) - 1][0]), "{:.2e}".format(y_n[len(y_n) - 1][1]), "{:.2e}".format(px), "{:.2e}".format(py)])
            
        linha = np.array(linha)
        print("m ----- n = 2^m ----- dt ----- \eta_x ----- \eta_y ---- px ----- py")
        print(linha)
