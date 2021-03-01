class MatrizDistancia:
    def __init__(self, arq):
        self.L = []                               # lista de nomes das OTUs:: L[i] = u <=> u eh a i-esima OTU em L
        self.H = {}                               # lista dos indices de L:: H[u] = i <=> L[i] = u
        self.n = 0                                # quantidade de OTUs  
        self.d = []                               # distancia entre OTUs:: matriz triangular. Para i > j, d[i][j] eh a distancia entre L[i] e L[j] 
        self.D = []                               # soma das distancias entre L[i] e as demais OTUs, i. e, toda a linha:: D[i] = \sum{j < i} d[i][j] + \sum{j > i} d[j][i]
        
        #### LEITURA DO ARQUIVO DE ENTRADA ####
        # ================================================================
        f = open(arq, "r")
        linha = f.readline().split()

        # definindo n
        self.n = len(linha)                       

        # definindo L e H
        for i in range(self.n):
            self.L.append(linha[i])
            self.H[linha[i]] = i

        # definindo a matriz triangular d
        for i in range(self.n): self.d.append([])

        i = 0
        while i < self.n:
            linha = f.readline().split()
            self.d[i] = list(map(float, linha[1:i+1]))
            
            i = i + 1
        f.close()
        # ================================================================

        # definindo D
        
        # ================================================================
        self.D = [0] * self.n                     # soma de toda a linha i 

        for i in range (self.n):
            for j in range(0,i):         self.D[i] = self.D[i] + self.d[i][j]
            for j in range(i+1, self.n): self.D[i] = self.D[i] + self.d[j][i]
        # ================================================================

    def vazia(self):
        return (self.n == 0) 
        
    def __str__(self):
        # ================================================================
        # valor de n
        saida = "n = " + str(self.n) + "\n\n"
        # ================================================================

        # ================================================================
        # valor de L, H e d
        saida = saida + "matriz d (completa) = \n"
        for i in range(self.n): saida = saida + "\t" + self.L[i] + " "
        saida = saida + "\n"

        for i in range(len(self.d)):
            saida = saida + self.L[i] + "\t"
            for j in range (len(self.d[i])): saida = saida + str(self.d[i][j]) + "\t"
            saida = saida + "0.0" + "\t"
            j = i + 1
            while j < len(self.d):
                saida = saida + str(self.d[j][i]) + "\t"
                j = j + 1
            saida = saida + "\n"
        # =================================================================
        # matriz d bruta
        saida = saida + "\nMatriz d Bruta (Estrutura de dados usada)\n"
        for i in range(self.n): saida = saida + "\t" + str(self.L[i])
        for i in range(len(self.d)):
            saida = saida + "\n" + str(self.L[i]) + "\t"
            for j in range (len(self.d[i])): saida = saida + str(self.d[i][j]) + "\t"
        # =================================================================
                
        # =================================================================
        # =================================================================
        # hash H
        saida = saida + "\n\nhash H: "
        saida = saida + str(self.H)
        # =================================================================

        # =================================================================
        # =================================================================
        # matriz D
        saida = saida + "\n\nmatriz D: "
        saida = saida + str(self.D)
        # =================================================================


        return saida               

    def neighborJoining (self, u, v):
        # se n == 2, remove u e v da matriz e devolve (u, v, d[L[u]][L[v]])
        # se n > 2, insere uma nova OTU uv, e remove u e v da matriz e devolve (u, v, uv, d[L[u]][L[uv]], d[L[uv]][L[v]])
        
        # =================================================================
        # =================================================================
        # caso especial n = 2, nao precisa fazer nada
        if self.n == 2:
            self.n = 0
            return u, v, "raiz", 0.5 * (int(100 * self.d[1][0])/100.), 0.5 * (int(100 * self.d[1][0])/100.) 
        # =================================================================
        # =================================================================

        L = self.insereOTU(u, v)

        self.removeOTU(u)
        self.removeOTU(v)
        return L

    def insereOTU (self, u, v):
        novo = u + v
        self.d.append([0] * self.n)
        self.L.append(novo)
        h = self.H[novo] = self.n
        
        i = self.H[u]
        j = self.H[v]
        
        if i < j: i, j = j, i

        self.d[h][i] = abs(0.5 * self.d[i][j] + 1.0 * (self.D[i] - self.D[j])/(2 * (self.n - 2)))
        self.d[h][j] = abs(self.d[i][j] - self.d[h][i])

        for k in range(j):      self.d[h][k] = abs(0.5 * (self.d[j][k] + self.d[i][k] - self.d[i][j]))
        for k in range(j+1, i): self.d[h][k] = abs(0.5 * (self.d[k][j] + self.d[i][k] - self.d[i][j]))
        for k in range(i+1, self.n): self.d[h][k] = abs(0.5 * (self.d[k][j] + self.d[k][i] - self.d[i][j]))

        for k in range(self.n): self.D[k] += self.d[h][k]
        self.D.append(0)
        for k in range(self.n): self.D[self.n] += self.d[h][k]
            
        self.n = self.n + 1
        return u, v, novo, int(100*self.d[h][i])/100., int(100*self.d[h][j])/100.


    def removeOTU(self, u):
        self.n = self.n - 1
        i = self.H[u]
        h = self.n

        for k in range(i):
            self.D[k] = self.D[k] - self.d[i][k]
            self.d[i][k] = self.d[h][k]

        self.D[i] = self.D[h] - self.d[h][i]

        for k in range(i+1, self.n):
            self.D[k] = self.D[k] - self.d[k][i] 
            self.d[k][i] = self.d[h][k]

        del self.H[self.L[i]]    
        self.L[i] = self.L[h]
        self.L.pop()
        self.H[self.L[i]] = i
        self.d.pop()
        self.D.pop()    

    
    
class Floresta:
    pass
    #     # =================================================================
    #     # nos de T ja calculados.
    #     saida = saida + "\nArvore T:\n" + str(self.T)
    #     # =================================================================


class MatrizQ:
    # matriz Q referente a uma MatrizDistancia dist
    def __init__(self, dist):
        self.L = [[]] * dist.n  # lista dos nomes das OTUs da MatrizDistancia dist
        self.Q = []             # matriz Q
        self.n = dist.n         # quantidade de linhas da matriz
        
        for i in range(dist.n):
            self.L[i] = dist.L[i]
            
            self.Q.append([])
            # calcular Q[i][j]
            for j in range(i): self.Q[i].append((dist.n - 2) * dist.d[i][j] - (dist.D[i] + dist.D[j]))


    def __str__(self):
        saida = "Matriz Q\n"
        for i in range(self.n): saida = saida + "\t" + str(self.L[i])
        for i in range(len(self.Q)):
            saida = saida + "\n" + str(self.L[i]) + "\t"
            for j in range (len(self.Q[i])): saida = saida + str(self.Q[i][j]) + "\t"
        return saida

    def minimo(self):
        u = 1
        v = 0
        for i in range (len(self.Q)):
            for j in range (len(self.Q[i])):
                if self.Q[i][j] < self.Q[u][v]: u, v = i, j
        
        return self.L[u], self.L[v]


    def minimoModificado(self, p):
        # esse metodo obtem uma lista R de \lceil self.n/p \rceil pares de OTUs que possuem menor valor Q
        
        # carregar a matriz Q em uma lista L = (i, j, Q[i][j]) onde i, j sao indices de Q
        L = []
        for i in range (len(self.Q)):
            for j in range (len(self.Q[i])):
                L.append([i, j, self.Q[i][j]])

        # transformar L em um heap        
        for h in range(1,len(L)):
            k = h
            while k > 0 and L[k][2] < L[(k-1)/2][2]:
                L[k], L[(k-1)/2] = L[(k-1)/2], L[k]
                k = (k-1)/2

        # indices eh uma lista dos indices das OTUs
        indices = []
        for i in range (len(self.Q)): indices.append(False)        
                
        # obter a lista de nos l[u] = True se o no u foi escolhido; False caso contrario.
        R = []
        i = 0
        while 2 * p * len(R) < len(self.Q):
            if indices[L[i][0]] == False and indices[L[i][1]] == False:
                R.append([self.L[L[i][0]], self.L[L[i][1]]])
                indices[L[i][0]] = indices[L[i][1]] = True
            i = i + 1
        
        # print "L = "
        # print L
        # print len(self.Q)
        print "R = "
        print R
        return R
 

T = [] # arvore
d = MatrizDistancia("matrizDistancia3.ent")
p = 2

while not d.vazia():
   Q = MatrizQ (d)
   R = Q.minimoModificado(p)

   while len(R) > 0:
       x, y = R.pop()
       T.append(d.neighborJoining(x, y))

   print Q
   print d
   print T

