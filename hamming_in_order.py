from queue import PriorityQueue
import numpy as np
from PIL import Image as im
import math

class HammingInOrderGraph:
    def __init__(self, r_val):
        #calculate variables based on r
        self.r = r_val
        self.n = pow(2, r_val) - 1
        self.k = self.n - self.r
        self.t = int((self.n * (self.n - 1)) / 6)
        self.P_count = pow(2, self.k)
        print(f'Hamming In-Order Graph (r = {self.r}):\t(n,k,d) = ({self.n},{self.k},3)\t|M| = {self.t}')
        
        #add placeholders (save work in case we don't need it)
        self.H = []
        self.Hp = []
        self.swaps = []
        self.G = []
        self.Gp = []
        self.M = []
        self.V = []
        self.E = []
    
    def buildH(self):
        #generate H (parity check matrix)
        self.H = []
        for i in range(1, self.n + 1):
            self.H.append(self.bin_str(i, self.r))
        
        #generate H' (parity check matrix in standard form) by performing row swaps
        self.Hp = []
        self.swaps = []
        for row in self.H:
            self.Hp.append(row)
        for i in range(self.r):
            start = pow(2, self.r - 1 - i) - 1
            #check value of ith row - have we already swapped?
            if int(self.Hp[i],2) != pow(2, self.r - 1 - i):
                from_row = start
                to_row = i
                while True:
                    self.swaps.insert(0,(to_row,from_row))
                    self.Hp[to_row] = self.H[from_row] #perform swap
                    from_row = to_row #increment from_row
                    if self.isPowerOfTwo((val := int(self.H[to_row],2))):
                        to_row = self.r - 1 - int(math.log2(val)) #increment to_row
                    else:
                        break
                self.Hp[start] = self.H[from_row] #perform final swap
                self.swaps.insert(0,(start,from_row))
    
    def printH(self):
        if self.H == []:
            self.buildH()
        self.printListAsMatrix(self.H,"H")
        self.printListAsMatrix(self.Hp,"H'")
    
    def buildG(self):
        if self.H == []:
            self.buildH()
        
        #generate G' (generating matrix in standard form)
        self.Gp = []
        for i in range(self.k):
            word = self.Hp[i + self.r] # this row of A
            for j in range(i):
                word += "0" # left padding for this row of I
            word += "1"
            while len(word) < self.n:
                word += "0" # right padding for this row of I
            self.Gp.append(word)
        
        #generate G (generating matrix) by performing column swaps
        self.G = []
        for row in self.Gp:
            self.G.append(row)
        for (from_col, to_col) in self.swaps:
            for i in range(self.k):
                self.G[i] = self.G[i][:to_col] + self.Gp[i][from_col] + self.G[i][(to_col + 1):]
    
    def printG(self):
        if self.G == []:
            self.buildG()
        self.printListAsMatrix(self.Gp,"G'")
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
        data.save('hamming_io_graph_' + str(self.r) + '.png')
    
    def weight(self, x):
        w = 0
        b = str(bin(x)) # binary string
        for c in b:
            if c == "1":
                w += 1
        return w
    
    def isPowerOfTwo(self, num):
        power = int (math.log(num, 2) + 0.5)
        return (2 ** power) == num
    
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
    hg = HammingInOrderGraph(r)
    hg.test()

# Test with r=3
testR(3)

# Test with r=4
testR(4)