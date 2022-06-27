import matplotlib.pyplot as plt
import os

def odczyt():
    with open("ksiazkauu.lis","rb") as f:
        data=f.read()
        f.close()
    return data

def wielkosc(x):
    i=0
    j=0
    w=[]
    while True:
        if(i<=3):
            w.append(int(x[0]))
            del x[0]
            i=i+1
        else:
             break
    return x,w
def pelna_skala(x):
    h=x[0]+x[1]*255
    w=x[2]+x[3]*255
    return h,w

def odbudowa_ksiazki_etap1(data):
    ii=True
    ile=0
    pamietaj=data[0]
    del data[0]
    x=[]
    y=[]
    k=""
    u=0
    for i in data:
        if ii:
            x.append(i)
            ii=False
            pamietaj=pamietaj-1
            u=u+2
        elif ii==False and ile==0:
            ile=i
            u=u+i
        else:
            n="{0:b}".format(i)
            
            for j in range (len(n)):
                if j==0:
                    continue
                else:
                    k=k+n[j]
            
            ile=ile-1
            if ile==0:
                y.append(k)
                k=''
                ii=True
            if(pamietaj==0):
                return x,y,u

def HuffmanCodeOdczyt(d):
    file=open("lisek.lis","wb")
    file.write(d)
    file.close()
    
    lista=[]
    file=open("lisek.lis","rb")
    znak=file.read(1)
    while znak:
        lista.append(int.from_bytes(znak,"little"))
        znak=file.read(1)
        
    x=""
    file.close()
    os.remove("lisek.lis")
    for i in range(0,len(lista)-1):
        konwersja=format(lista[i],"b").zfill(8)
        x=x+konwersja
    
    #x=x+format(lista[len(lista)-2],"b").zfill(8-lista[len(lista)-1])
    #if len(x)%8!=0:
     # for i in range (8-len(x)%8):
      #  x=x+"0"
    return x
    
#TODO - HuffmanDecode
def HuffmanDecode(CodeBook, Code):
  dekodowanie=[]
  x=""
  for i in range (len(Code)):
    x=x+Code[i]
    for j in CodeBook.keys():
        if CodeBook[j]==x:
            x=''
            dekodowanie.append(int(j))
            break

    
  return dekodowanie    

def ksiazka_z_pliku_reaktywacja(znak,kod):
    ksiega={}
    for i in range (len(znak)):
        s=str(znak[i])
        ksiega[s]=kod[i]
    return ksiega

def zostaw_w_pliku_tylko_obraz(u,x):
    for i in range (u):
        del x[0]
    return x

def dispImagesWithHistograms(image1,image2):
  fig,axs = plt.subplots(2,2, figsize=(12,12))
  axs[0][0].imshow(image1, cmap=plt.cm.gray,vmin=0, vmax=255)
  axs[1][0].hist(image1.ravel(),bins=256,range=[0,255])
  axs[0][1].imshow(image2, cmap=plt.cm.gray,vmin=0, vmax=255)
  axs[1][1].hist(image2.ravel(),bins=256,range=[0,255])
  fig.tight_layout()
  return