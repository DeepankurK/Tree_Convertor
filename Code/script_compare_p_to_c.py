import os
for p in range(1,11):
    file= open(r"/home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/noisy_genotype_dataset"+str(p)+".txt","r")
    content=file.readlines()
    for j,i in enumerate(content):
        a=i.split(" ")
        del a[0]
        for k in range(0,len(a)):
            if a[k]=='2':
                a[k]='1'
        content[j]=" ".join(a)
    file.close()
    file=open("/home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/dataset"+str(p)+".txt",'w')
    file.write("".join(content))
    file.close()
    file=open("/home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/Orig_tree_dataset"+str(p)+".txt",'r')
    content=file.readlines()
    file.close()
    #if len(content)>1:
     #   content=content[-1]
    #print(content[1],len(content))
    
    file=open("/home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/orig"+str(p)+".txt",'w')
    file.write(content[1])
    file.close()
    for j in range(1,7):
        os.system("echo algo "+str(j)+" >> result"+str(p)+".txt")
        for k in [1,2]:
            os.system("python compare.py -t 2 -u 1 -a /home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/orig"+str(p)+".txt -b /home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/noisy_genotype_dataset"+str(p)+"_mlTree.newick -m "+str(k)+" -x -s /home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/test_Data/10_clones/dataset"+str(p)+"/dataset"+str(p)+".txt -r 2 -l "+str(j)+" -o  /home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/results_p_to_c/"+str(p)+"_"+str(j)+"_"+str(k)+" >> result"+str(p)+".txt")
            print(p,j,k)
            
import numpy as np
ls=np.zeros((10,12))
for p in range(1,11):
    file=open("/home/deepank/Downloads/Prof_Hamim/Convertor_project/Code/result"+str(p)+".txt",'r')
    content=file.readlines()
    k=0
    for i in content:
        print(k,i)
        if 'core' in i:
            if 'pairwise' in i:
                #print(i.split(':')[1].strip())
                ls[p-1][2*k]=float(i.split(':')[1].strip())
            elif 'Rand Score' in i:
                #print(i.split(':')[1].strip(),i)
                ls[p-1][2*k+1]=float(i.split(':')[1].strip())
                k+=1
    file.close()
print(ls)
import pandas as pd
df=pd.DataFrame(ls)
print(df)    
df.to_csv('result_p_to_c.txt',index=False,header=False)
