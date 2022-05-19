
import pandas as pd

file= open(r'/home/deepank/Downloads/Prof_Hamim/Convertor_project/dataset2.txt')
content=file.readlines()

for j,i in enumerate(content):
    a=i.split(" ")
    a[len(a)-1]=a[len(a)-1].split('\n')[0]
    content[j]=" ".join(a)
df=pd.DataFrame(content,columns=['data'])
df.to_csv('scite_dataset2.csv',index=False,header=False)
