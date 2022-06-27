import numpy as np
import math 
import collections
import heapq


#koder
def DeQuant2(Value, Nbits, Range) :
  W = 2**(Nbits-1)
  T = Range*(2**(1-Nbits))
  return (Value+W)*T - Range
#---------------------------------
def Quant2(Value, Nbits, Range) :
  SignMax = 2**(Nbits-1)-1
  SignMin = -(1+SignMax)
  Q   = 2**Nbits-1
  T   = (Q+1)/(2*Range)
  u   = (Value+Range)*T
  u   = u+SignMin
  #u(u < SignMin) = SignMin
  #u(u > SignMax) = SignMax 
  return math.floor(u)
#--------------------------------
def NonUniQuant2(Value, Nbits, Range, Mu) :
  Value = np.double(Value)
  NUQ   = Range * (math.log10(1+Mu*np.abs(Value)/Range)/math.log10(1+Mu))*np.sign(Value)
  return np.int16(Quant2(NUQ, Nbits, Range))
#--------------------------------
def NonUniDeQuant2(Value, Nbits, Range, Mu) :
  Value = np.double(Value)
  NUDQ = (Range / Mu) * (10**((math.log10(1+Mu)/Range)*np.abs(Value))-1)*np.sign(Value);
  return np.int16(DeQuant2(NUDQ, Nbits, Range))
#---------------------------------
def Dequant(error_q, pred, nbits, qnu) : 
  if nbits == 0 and qnu == 0 :   
    error_deq = error_q
    recon     = pred + error_deq
  elif nbits != 0 and qnu == 0 :   
    error_deq = DeQuant2(error_q,nbits,255)
    recon = pred + error_deq
  else :
    error_deq = NonUniDeQuant2(error_q,nbits,255,qnu)
    recon = pred + error_deq
  return recon
#---------------------------------
def QuantDequant(org, pred, nbits, qnu) :
  if nbits == 0 and qnu == 0 :
    error     = org - pred
    error_q   = error
    error_deq = error_q
    recon     = pred + error_deq
  elif nbits !=0 and qnu == 0 :
    error     = org - pred
    error_q   = Quant2(error,nbits,255)
    error_deq = DeQuant2(error_q,nbits,255)
    recon     = pred + error_deq
  else :
    error     = org - pred
    error_q   = NonUniQuant2(error,nbits,255,qnu)
    error_deq = NonUniDeQuant2(error_q,nbits,255,qnu)
    recon     = pred + error_deq
  return error_q, recon
def Predcode (img_org, pred_type, nbits=0, qnu=0) : 
  img_org_np = np.array(img_org,dtype=np.int16)
  dim          = img_org_np.shape
  error_mat    = np.zeros_like(img_org_np)
  img_rec      = np.zeros_like(img_org_np)


  if pred_type == 1 : 
    for y in range(0,dim[0]) :      #loop over 
      if y == 0 :                     # first line - first point
        pred                           = 128
        error_mat[y,0], img_rec[y,0]   = QuantDequant(img_org_np[y,0], pred, nbits, qnu)
      else :                        # remaining lines - first point
        pred = img_rec[y-1,0]
        error_mat[y,0], img_rec[y,0]   = QuantDequant(img_org_np[y,0], pred, nbits, qnu)
      for x in range(1,dim[1]) :    # loop over collumns
        pred                           = img_rec[y,x-1]
        error_mat[y,x], img_rec[y,x]   = QuantDequant(img_org_np[y,x], pred, nbits, qnu)
  elif pred_type == 3 :
    for y in range(0,dim[0]) :      #loop over lines
      if y == 0 :                     # first line - first point
        pred                           = 128
        error_mat[y,0], img_rec[y,0]   = QuantDequant(img_org_np[y,0], pred, nbits, qnu)
        for x in range(1,dim[1]) :  # loop over columns
          pred                         = img_rec[y,x-1]
          error_mat[y,x], img_rec[y,x] = QuantDequant(img_org_np[y,x], pred, nbits, qnu)
      else :                        # remaining lines - first point
        pred                           = img_rec[y-1,0]
        error_mat[y,0], img_rec[y,0]   = QuantDequant(img_org_np[y,0], pred, nbits, qnu)
        for x in range(1,dim[1]) :  #loop over collumns
          pred                         = np.int16(np.double(img_rec[y,x-1] + img_rec[y-1,x-1] + img_rec[y-1,x])/3.0)
          error_mat[y,x], img_rec[y,x] = QuantDequant(img_org_np[y,x], pred, nbits, qnu)
  else :                  # pred_type not recognized
    print('\nOnly left (1) and three-points (3) prediction is supported!\n\n')
  
  if nbits == 0 and qnu == 0 :
    error_img = np.uint8(error_mat+128)
  elif nbits != 0 and qnu == 0 :
    error_img = np.uint8(DeQuant2(error_mat,nbits,255)+128);
  else :
    error_img = np.uint8(NonUniDeQuant2(error_mat,nbits,255,qnu)+128)

  return error_img
#-----------------------------------------------------------------------------------------------------------------------
# Huffmann codig functions
#-----------------------------------------------------------------------------------------------------------------------
class xHuffTree:
  def __init__(self, Value, Count, Left=None, Right=None):
    if(Value != None):
      self.Value = Value
      self.Count = int(Count)
      self.Left  = None
      self.Right = None
    else:
      self.Value = None
      self.Count = Left.Count + Right.Count
      self.Left  = Left
      self.Right = Right
    return  
  def __lt__(self, other): return (self.Count < other.Count)

def __init__(self):
  self.HuffmanTree = []
  return

#-----------------------------------------------------------------------------------------------------------------------

    
def xHuffmanTree(Vals):
  HuffmanTree = []
  #insert values
  for Val, Count in Vals.items():
    heapq.heappush(HuffmanTree, xHuffTree(Val, Count, None, None))
  #build tree
  while len(HuffmanTree)>1:
    L = heapq.heappop(HuffmanTree)
    R = heapq.heappop(HuffmanTree)
    heapq.heappush(HuffmanTree, xHuffTree(None, None, L,R))
  return HuffmanTree[0]

#-----------------------------------------------------------------------------------------------------------------------

def xGenerateCode(Node:xHuffTree, Code, Length):
  CodeTable = {}
  if(not Node.Left and not Node.Right):
    CodeTable[Node.Value] = Code
  else:
    CodeTable.update(xGenerateCode(Node.Left , Code+"0", Length+1))
    CodeTable.update(xGenerateCode(Node.Right, Code+"1", Length+1))
  return CodeTable



#-----------------------------------------------------------------------------------------------------------------------

def HuffmanCodeBook(Vals):
  
  Tree     = xHuffmanTree(Vals)
  CodeBook = xGenerateCode(Tree, "", 0)
  return CodeBook

#-----------------------------------------------------------------------------------------------------------------------

def HuffmanEncode(CodeBook, Data):
  CodeSymbols = [CodeBook[D] for D in Data]
  Code = "".join(CodeSymbols)
  return Code

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

def HuffmanCodeBits(ImgData):
  DataIn  = np.array(ImgData).ravel()
  Counts  = collections.Counter(DataIn)
  Book    = HuffmanCodeBook(Counts)
  Code    = HuffmanEncode(Book, DataIn)
  CodeLen = len(Code)
  DataLen = len(DataIn)
  BitsPerSample = CodeLen / DataLen
  return BitsPerSample

#-----------------------------------------------------------------------------------------------------------------------

def HuffmanCodeZapis(kod):
    x=bytearray()
    for i in range(0,len(kod),8):
        x.append(int(kod[i:i+8],2))
    ile= 8- len(kod)%8
    x.append(ile)
    return x


def odczyt_ksiazki_z_pliku():
  with open("ksiazka.lis","rb") as f:
    data=f.read()
  return data
  
def pred(img):
    tryb_predyktora1_3 = 3
    predykcja =Predcode(img,tryb_predyktora1_3,8) 
    return predykcja
    
def symbole_ksiazka(predykcja):
    x=[]
    y=[]
    x.append(predykcja[0][0])
    y.append(0)
    t=False
    for i in range (len(predykcja)):
        for j in range (512):
            for k in range(len(x)):
                if x[k]==predykcja[i][j]:
                    y[k]=y[k]+1
                    t=True
            if t==False:
                x.append(predykcja[i][j])
                y.append(1)
            t=False
    return x,y
    
def liczba_symboli_calkowita(y):
    k=0        
    for i in range(len(y)):
        k=k+y[i]
    return k
    
def kodowanie(p,book):
    ciag=""
    for i in range(len(p)):
        temp=str(p[i])
        ciag=ciag+book[temp]
    return ciag

def ksiazka_prawdopodobienst(x,y):
    symb={}
    k=liczba_symboli_calkowita(y)
    for i in range (len(x)):
        tem=str(x[i])
        symb[tem]=y[i]/k
    return symb
def na_bin(kk):
    kk='1'+kk
    x=len(kk)
    #print(kk)
    k= x-1
    w=0
    for i in range (x):
        w=w+int(kk[i])*2**k
        k=k-1


    return w 


def ksiazka_bit(ksiazka):
    liczba=0
    znaki=[]
    temp=[]
    bity=[]
    for i in ksiazka.keys():
        x=""
        znaki.append(i)
        k=ksiazka[i]
        for j in range(len(ksiazka[i])):
            x=x+k[j]
            if len(x)==6:
                liczba=liczba+1
                bity.append(na_bin(x))
                x=''
        if len(x)>0:
            liczba=liczba+1
            bity.append(na_bin(x))
            x=''
        temp.append(liczba)
        liczba=0

    return znaki,temp,bity 
def podziel(x):
    pamietaj=0
    while True:
        if x>255:
            x=x-255
            pamietaj=pamietaj+1
        else:
            break
    return x,pamietaj
    
def dopisz_ksiazke(znaki,temp,bity):
    b=bytearray()
    cale=0
    for i in range(len(znaki)):
        z=znaki[i]
        b.append(int(z))
        z=temp[i]
        b.append(int(z))
        for j in range(z):
            c=bity[cale]
            cale=cale+1
            b.append(int(c))
    return b
def zapis_pliku(k):
    file=open("ksiazkauu.lis","wb")
    file.write(k)
    file.close()
def zapis_do_pliku_rozmiaru(x1,y1,x2,y2):
    b=bytearray()
    b.append(int(x1))
    b.append(int(y1))
    b.append(int(x2))
    b.append(int(y2))

    return b
