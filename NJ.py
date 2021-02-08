class MatrizDistancia:
    def __init__(self, arq):
        self.T = []                               # arvore construida. Lista: cada elemento possui x, y, c
                                                  # onde x e y sao nos vizinhos em T e c eh a distancia entre x e y 
        
        self.L = []                               # lista de nomes dos AA
        self.H = {}                               # lista dos indices de L
        
        f = open(arq, "r")
        linha = f.readline().split()
        for i in range(len(linha)):
            self.L.append(linha[i])
            self.H[linha[i]] = i

        self.n = len(linha)                        # quantidade de sequencias #
        # print self.n
        
        # ================================================================
        # ================================================================
        self.d = []                               # distancia entre os AA
        for i in range(self.n): self.d.append([])

        # note que d vai armazenar somente uma matriz triangular
        i = 0
        while i < self.n:
            linha = f.readline().split()
            self.d[i] = list(map(float, linha[1:i+1]))
            
            i = i + 1
        f.close()
        # ================================================================
        # ================================================================
        self.D = [0] * self.n                     # soma de toda a linha i 
        i = 0
        while i < self.n:
            # print i
            j = 0
            while j < i:
                self.D[i] = self.D[i] + self.d[i][j]
                j = j + 1
            j = i + 1
            while j < self.n:
                self.D[i] = self.D[i] + self.d[j][i]
                j = j + 1
            i = i + 1
        # ================================================================
            
        
        
    def __str__(self):
        # ================================================================
        # valor de n
        saida = "n = " + str(self.n) + "\n\n"
        # ================================================================

        # ================================================================
        # valor de L e d
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
        # matriz D
        saida = saida + "\n\nmatriz D: "
        saida = saida + str(self.H)
        # =================================================================
                
        # =================================================================
        # nos de T ja calculados.
        saida = saida + "\nArvore T:\n" + str(self.T)
        # =================================================================

        return saida               

    def insereOTU (self, u, v): # insere uma unica OTU
        # caso especial n = 2, nao precisa fazer nada
        if self.n == 2: 
            lista = u, v, int(100 * self.d[1][0])/100.
            self.T.append(lista)
            self.n = 0
            return self.d[1][0]
        # =================================================================
                
        # =================================================================
        # caso geral
        novo = u + v
        self.d.append([0] * self.n)
        self.L.append(novo)
        h = self.H[novo] = self.n
        print u, v, novo
        
        i = self.H[u]
        j = self.H[v]
        
        if i < j: i, j = j, i

        self.d[h][i] = 0.5 * self.d[i][j] + 1.0 * (self.D[i] - self.D[j])/(2 * (self.n - 2))
        self.D[i] += self.d[h][i]
        # lista = u, novo, self.d[h][i]
        # self.T.append(lista)
        
        self.d[h][j] = self.d[i][j] - self.d[h][i]
        self.D[j] += self.d[h][j]
        # lista = v, novo, self.d[h][j]
        lista = u, v, novo, int(100*self.d[h][i])/100., int(100*self.d[h][j])/100.
        self.T.append(lista)
                
        k = 0
        while k < j:
            self.d[h][k] = 0.5 * (self.d[j][k] + self.d[i][k] - self.d[i][j])
            self.D[k] += self.d[h][k]
            k = k + 1
        k = k + 1
        while k < i:
            self.d[h][k] = 0.5 * (self.d[k][j] + self.d[i][k] - self.d[i][j])
            if self.d[h][k] < 0: self.d[h][k] *= -1
            self.D[k] += self.d[h][k]
            k = k + 1
        k = k + 1
        while k < self.n:
            self.d[h][k] = 0.5 * (self.d[k][j] + self.d[k][i] - self.d[i][j])
            self.D[k] += self.d[h][k]
            k = k + 1
            
        self.n = self.n + 1

        self.D.append(0)
        for k in range(self.n - 1): self.D[self.n - 1] += self.d[h][k]

        return self.d[i][j] 

    def removeOTU(self, u):
        if self.n == 0: return
        self.n = self.n - 1
        i = self.H[u]
        h = self.H[self.L[self.n]]

        j = 0
        while j < i:
            self.D[j] = self.D[j] - self.d[i][j]
            self.d[i][j] = self.d[h][j]
            j = j + 1
        print i, h, self.n    
        self.D[i] = self.D[h] - self.d[h][i]
        j = j + 1
        while j < self.n:
            self.D[j] = self.D[j] - self.d[j][i] 
            self.d[j][i] = self.d[h][j]
            j = j + 1

        del self.H[self.L[i]]    
        self.L[i] = self.L[h]
        self.L.pop()
        self.H[self.L[i]] = i
        self.d.pop()
        self.D.pop()

    
        
class MatrizQ:
    def __init__(self, dist):
        self.Q = []
        self.L = [[]] * dist.n
        for i in range(dist.n):
            self.L[i] = dist.L[i]
            self.Q.append([])
            for j in range(i):
                # calcular Q[i][j]
                self.Q[i].append((dist.n - 2) * dist.d[i][j] - (dist.D[i] + dist.D[j]))
        # print "Q = "
        # for i in range(dist.n): print self.Q[i]

    def minimoQ(self):
        u = 1
        v = 0
        for i in range (len(self.Q)):
            for j in range (len(self.Q[i])):
                if self.Q[i][j] < self.Q[u][v]: u, v = i, j
        # print u, v
        return self.L[u], self.L[v]

    def minimoModificado(self, p):
        # esse metodo obter a lista R de \lceil self.n/p \rceil pares que possuem menor valor Q
        
        # carregar a matriz Q em uma lista
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

        indices = []
        for i in range (len(self.Q)): indices.append(False)        
                
        # obter a lista de nos l[u] = True se o no u foi escolhido; False caso contrario.
        R = []
        i = 0
        while len(R) < (len(self.Q) + p - 1)/p:
            if indices[L[i][0]] == False and indices[L[i][1]] == False:
                R.append([self.L[L[i][0]], self.L[L[i][1]]])
                indices[L[i][0]] = indices[L[i][1]] = True
            i = i + 1
        
        # print "L = "
        # print L
        # print len(self.Q)
        # print "R = "
        # print R
        return R
                       
        
    
d = MatrizDistancia("matrizDistancia3.ent")
soma = 0
print d
while d.n > 1:
    q = MatrizQ(d)
    R = q.minimoModificado(3)
    for i in range(len(R)):
        x, y = R[i][0], R[i][1]
        soma = soma + d.insereOTU (x, y)
    print d
    for i in range(len(R)):
        x, y = R[i][0], R[i][1]
        print "remove", x
        d.removeOTU (x)
        print "remove", y
        d.removeOTU (y)
    print d
print "Valor total: ", soma
