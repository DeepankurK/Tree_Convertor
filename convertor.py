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
import sys
import getopt

path_input=""
path_matrix=""
path_output=""
filename="Converted_Tree"
typ=0
dot=""
content=""

def graph_to_dot(content):
    c=0
    dot=Digraph()
    dot.attr('node')
    #print(dot)
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
                    if j+2<len(i):
                        dot.edge(str(int(i[j-3]+i[j-2])),str(int(i[j+1]+i[j+2])))
                    else:
                        dot.edge(str(int(i[j-3]+i[j-2])),str(int(i[j+1])))
    #print(dot
    return dot  
  
def dot_to_newick(dot):
    label={}
    #print(dot)
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

def input_func(ch,content,gr=False):
    label={}
    if ch=='c' and gr:       
        dot=graph_to_dot(content)
        content,label=dot_to_newick(dot)
    handle=io.StringIO(content)
    tree=Phylo.read(handle,"newick")
    Phylo.draw(tree,do_show=False)
    plt.savefig(path_output+'Original_tree.png')
    if ch=='c':
        return tree,matrix,label
    return tree,matrix
def write(content):
    file=open(path_output+'converted_newick.txt','w')
    file.write(content)
    file.close()
    
def phylo_tree(content,typ=0):
    in_tree,matrix=input_func('p',content)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 2 for Clonal Tree \n 3 for Mutation Tree \n"))
    if typ==3:
        obj=phylo.Phylo_to_Mut(in_tree,matrix,z)
    if typ==2:
        obj=phylo.Phylo_to_Clonal(in_tree,matrix,z)
    dot=obj.convert()
    #print(dot)
    newick,_=dot_to_newick(dot)
    return newick,dot
    
def clonal_tree(gr,content,typ=0):
    in_tree,matrix,label=input_func('c',content,gr)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 1 for Phylogenetc Tree \n 3 for Mutation Tree \n"))
    if typ==1:
        obj=clonal.Clonal_to_Phylo(in_tree,matrix,label,gr,z)
    if typ==3:
        obj=clonal.Clonal_to_Mut(in_tree,matrix,label,gr,z)
    dot=obj.convert()
    newick,_=dot_to_newick(dot)
    return newick,dot
    
def muta_tree(content,typ=0):
    in_tree,matrix=input_func('m',content)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 1 for Phylogenetic Tree \n 2 for Clonal Tree \n"))
    if typ==2:
        obj=mut.Mut_to_Clonal(in_tree,matrix,z)
    if typ==1:
        obj=mut.Mut_to_Phylo(in_tree,matrix,z)
    dot=obj.convert()
    newick,_=dot_to_newick(dot)
    return newick,dot

gr=False
z=False
arg=sys.argv
if len(arg)==1:
    ch=int(input("Enter the tree format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\n"))
    gr=bool(input("The input tree is in igraph format? 1 or 0?\n"))
    x=bool(input("The numbering of cell/mutations in the current tree starts from 0 or 1?\n"))
    
else:
    arg=arg[1:]
    short_options="t:i:m:o:f:r:gz"
    long_options=["tree","input","matrix","output","filename","result","igraph","zero"]
    try:
        arguments, values = getopt.getopt(arg, short_options, long_options)
    except getopt.error as err:
        print (str(err))
        sys.exit(2)
    for curr_arg, curr_val in arguments:
        if curr_arg in ("-t","--tree"):
            ch=int(curr_val)
        if curr_arg in ("-i","--input"):
            path_input=curr_val
        if curr_arg in ("-m","--matrix"):
            path_matrix=curr_val
        if curr_arg in ("-o","--output"):
            path_output=curr_val
        if curr_arg in ("-f","--filename"):
            filename=curr_val
        if curr_arg in ("-r","--result"):
            typ=int(curr_val)
        if curr_arg in ("-g","--igraph"):
            gr=True
        if curr_arg in ("-z","--zero"):
            z=True
#print(gr)
if path_output=="":
    if not os.path.exists(os.getcwd()+'/results/'):
        os.makedirs(os.getcwd()+'/results/')
    path_output=os.getcwd()+'/results/'
else:
    if path_output[-1]!='/':
        path_output+='/'
    if not os.path.exists(path_output+'results/'):
        os.makedirs(path_output+'/results/')
    path_output+='results/'
if path_input=="":
    path_input=input("Enter the path of input file.\n")
content=open(path_input,"r").read()
if path_matrix=="":
    path_matrix=input("Enter the path of cell mutation matrix file.\n")
matrix=pd.read_csv(path_matrix,header=None,sep=" ")
matrix=matrix.replace(2,1)

print("Output will be stored in results folder.")
if ch==1:
    newick,dot=phylo_tree(content,typ)
elif ch==2:
    newick,dot=clonal_tree(gr,content,typ)
elif ch==3:
    newick,dot=muta_tree(content,typ)
else:
    print("Wrong Tree Choice inputted.")
    sys.exit(2)
dot.render(filename=filename,directory=path_output)
write(newick)
