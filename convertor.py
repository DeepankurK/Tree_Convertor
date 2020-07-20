from  newick import loads
from Bio import Phylo
import io
import pandas as pd
import re
import phylo
import mut
import clonal
import os
from ete3 import Tree
from graphviz import Digraph
import matplotlib.pyplot as plt
import shutil
#    cell file
 #   cell_file=input("Enter the path of file containing cell names or enter N for 1 to n to be considered.\n")

#output file
#terminal_getops
 
    
def graph_to_dot(content):
    c=0
    dot=Digraph()
    dot.attr('node')
    for i in content.split('\n'):
        if c!=0:
            q=i.strip().split(" ")
            dot.node(str(c)," ".join(q[1:]))
            c=0
        if 'clones' in i:
            c=re.sub('\D', '',str(i))
            c= int(c) if c!='' else 0
    for i in content.split('\n'):
        if "->" in i:
            for j,k in enumerate(i):
                if k=='>':
                    dot.edge(str(int(i[j-3]+i[j-2])),str(int(i[j+1]+i[j+2])))
    #print(dot
    return dot  
  
def dot_to_newick(dot):
    label={}
    print(dot)
    for i in str(dot).split("\n"):
        if "label=" in i:
            q=i.split('[label=')
            k=q[1].split('"')
            if len(k)==1:
                label[q[0].strip()]=k[0].split("]")[0]
            else:
                label[q[0].strip()]=k[1].strip()
    #print(label)
    #print("done")
    table=[]
    for i in str(dot).split("\n"):
        if "->" in i:
            q=i.split('->')
            table.append((q[0].strip(),q[1].strip(),0.5))
    #print(table)
    tree = Tree.from_parent_child_table(table)
    q=tree.write(format=1)
    for i in label.keys():
        if ','+str(i)+':' not in q and '('+str(i)+':' not in q and ')' +str(i)+':' not in q:
            q=q.split(';')[0]+str(i)+':0.5;'
    #print(q)
    return q,label

def input_func(ch):
    input_file=input("Enter the path of input file.\n")
    content=open(input_file,"r").read()
    if ch=='c':       
        dot=graph_to_dot(content)
        content,label=dot_to_newick(dot)
    in_tree=loads(content)
    matrix_file=input("Enter the path of cell mutation matrix file.\n")
    matrix=pd.read_csv(matrix_file,header=None,sep=" ")
    matrix=matrix.drop([0],axis=1)
    matrix=matrix.replace(2,1)
    print("Output will be stored in results folder.")
    if os.path.exists(os.getcwd()+'/results/'):
        shutil.rmtree(os.getcwd()+'/results/')
    os.makedirs(os.getcwd()+'/results/')
    
    handle=io.StringIO(content)
    tree=Phylo.read(handle,"newick")
    Phylo.draw(tree,do_show=False)
    plt.savefig(os.getcwd()+'/results/Original_tree.png')
    
    if ch=='c':
        return in_tree,matrix,label
    return in_tree,matrix

def write(content):
    file=open(os.getcwd()+'/results/newick.txt','w')
    file.write(content)
    file.close()
    
def phylo_tree():
    in_tree,matrix=input_func('p')
    ch=int(input("Select the tree for conversion:-  1 . Mutation Tree and 2. Clonal Tree.\n"))
    if ch==1:
        obj=phylo.Phylo_to_Mut(in_tree,matrix)
    if ch==2:
        obj=phylo.Phlyo_to_Clonal(in_tree,matrix)
    dot=obj.convert()
    newick,_=dot_to_newick(dot)
    dot.render(filename="Converted",directory=os.getcwd()+'/results/')
    write(newick)
    
def clonal_tree():
    in_tree,matrix,label=input_func('c')
    ch=int(input("Select the tree for conversion:-  1 . Phylogenetc Tree and 2. Mutation Tree.\n"))
    if ch==1:
        obj=clonal.Clonal_to_Phylo(in_tree,matrix)
    if ch==2:
        obj=clonal.Clonal_to_Mut(in_tree,matrix)
    dot=obj.convert()
    newick,_=dot_to_newick(dot)
    dot.render(filename="Converted",directory=os.getcwd()+'/results/')
    write(newick)
    
def muta_tree():
    in_tree,matrix=input_func('m')
    ch=int(input("Select the tree for conversion:-  1 . Clonal Tree and 2. Phylogenetic Tree.\n"))
    if ch==1:
        obj=mut.Mut_to_Clonal(in_tree,matrix)
    if ch==2:
        obj=phylo.Mut_to_Phylo(in_tree,matrix)
    dot=obj.convert()
    newick,_=dot_to_newick(dot)
    dot.render(filename="Converte",directory=os.getcwd()+'/results/')
    write(newick)

ch=int(input("Enter the tree format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\n"))
if ch==1:
    phylo_tree()
elif ch==2:
    clonal_tree()
elif ch==3:
    muta_tree()
else:
    print("Wrong Choice inputted.")
    
    

