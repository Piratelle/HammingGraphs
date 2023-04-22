from queue import PriorityQueue
import numpy as np
from PIL import Image as im
import math

class HammingGraph:
    def __init__(self, r_val):
        #calculate variables based on r
        self.r = r_val
        self.n = pow(2, r_val) - 1
        self.k = self.n - self.r
        self.t = int((self.n * (self.n - 1)) / 6)
        self.P_count = pow(2, self.k)
        print(f'Hamming Graph (r = {self.r}):\t(n,k,d) = ({self.n},{self.k},3)\t|M| = {self.t}')
        
        #add placeholders (save work in case we don't need it)
        self.H = []
        self.G = []
        self.M = []
        self.V = []
        self.E = []
    
    def buildH(self):
        #generate words for H
        hq = PriorityQueue()
        for i in range(1, self.n + 1):
            w = self.weight(i)
            priority = ((w - 1) * self. n) - i
            hq.put((priority, self.bin_str(i, self.r)))
        
        #generate H (parity check matrix)
        self.H = []
        while not hq.empty():
            self.H.append(hq.get()[1])
    
    def printH(self):
        if self.H == []:
            self.buildH()
        self.printListAsMatrix(self.H,"H")
    
    def buildG(self):
        if self.H == []:
            self.buildH()
        
        #generate G (generating matrix)
        self.G = []
        for i in range(self.k):
            word = self.H[i + self.r] # this row of A
            for j in range(i):
                word += "0" # left padding for this row of I
            word += "1"
            while len(word) < self.n:
                word += "0" # right padding for this row of I
            self.G.append(word)
    
    def printG(self):
        if self.G == []:
            self.buildG()
        self.printListAsMatrix(self.G,"G")
    
    def buildM(self):
        if self.H == []:
            self.buildH()
        
        #generate M (set of non-zero min weight codewords)
        self.M = []
        for r1 in range(self.n):
            for r2 in range(r1 + 1, self.n):
                x1 = int(self.H[r1],2)
                x2 = int(self.H[r2],2)
                r3 = self.H.index(self.bin_str(x1 ^ x2, self.r))
                if r3 > r2:
                    v = [0 for i in range(self.n)]
                    v[r1] = 1
                    v[r2] = 1
                    v[r3] = 1
                    self.M.insert(0,''.join(map(str,v)))
    
    def printM(self):
        if self.M == []:
            self.buildM()
        self.printListAsMatrix(self.M,"M")
    
    def buildV(self):
        if self.G == []:
            self.buildG()
        
        #generate V (vertex set)
        self.V = []
        for i in range(self.P_count):
            w = self.bin_str(i, self.k)
            x = 0
            for j in range(self.k):
                if w[j] == "1":
                    x = x ^ int(self.G[j],2)
            self.V.append(x)
    
    def printV(self):
        if self.V == []:
            self.buildV()
        
        print() # print blank line for spacing
        for i in range(len(self.V)):
            print(f'v{i} =\t{self.bin_str(i, self.k)} * G =\t{self.bin_str(self.V[i], self.n)}')
    
    def buildE(self):
        if self.V == []:
            self.buildV()
        
        #generate E (edge set) as adjacency matrix
        self.E = [[0 for j in range(self.P_count)] for i in range(self.P_count)]
        for i in range(self.P_count):
            for j in range(i + 1, self.P_count):
                u = self.V[i]
                v = self.V[j]
                if self.weight(u ^ v) == 3:
                    self.E[i][j] = 1
                    self.E[j][i] = 1
    
    def printE_Mtx(self):
        if self.E == []:
            self.buildE()
        edges = []
        for i in range(self.P_count):
            edges.append(''.join(map(str,self.E[i])))
        self.printListAsMatrix(edges, "Edges")
    
    def printE(self):
        if self.E == []:
            self.buildE()
        print() # print blank line for spacing
        #print as neighbor set
        for i in range(self.P_count):
            neighbors = ""
            for j in range(self.P_count):
                if self.E[i][j] == 1:
                    neighbors += " v" + str(j)
            print(f'v{i}: ({neighbors} )')
    
    def drawE(self):
        if self.E == []:
            self.buildE()
        array = np.array(self.E).astype(np.uint8)
        for i in range(self.P_count):
            for j in range(self.P_count):
                array[i][j] *= 255
        #print(array)
        data = im.fromarray(array,'L')
        data.save('hamming_graph_' + str(self.r) + '.png')
    
    def weight(self, x):
        w = 0
        b = str(bin(x)) # binary string
        for c in b:
            if c == "1":
                w += 1
        return w
    
    def bin_str(self, x, length):
        s = str(bin(x))
        s = s[2:]
        while len(s) < length:
            s = "0" + s
        return s
    
    def printListAsMatrix(self, lst, lbl):
        print() # print blank line for spacing
        r = len(lst)
        for i in range(r):
            word = ""
            for c in lst[i]:
                word += " " + c
            if i < 1:
                print(f'{lbl}:\t⌈{word} ⌉')
            elif i < r - 1:
                print(f'\t|{word} |')
            else:
                print(f'\t⌊{word} ⌋')
    
    def test(self):
        self.printH()
        self.printG()
        self.printM()
        self.printV()
        if self.r < 4:
            self.printE_Mtx()
        self.printE()
        self.drawE()

def testR(r):
    hg = HammingGraph(r)
    hg.test()

# Test with r=3
testR(3)

# Test with r=4
testR(4)