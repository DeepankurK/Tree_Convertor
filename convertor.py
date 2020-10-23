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

path_output=""
path_matrix=""
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
  
def dot_to_newick(dot,save=False):
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
            if 's' not in q[1].strip().split(';')[0]:
                table.append((q[0].strip(),q[1].strip().split(';')[0],0.5))
    #print(table)
    tree = Tree.from_parent_child_table(table)
    q=tree.write(format=1)
    for i in label.keys():
        if ','+str(i)+':' not in q and '('+str(i)+':' not in q and ')' +str(i)+':' not in q:
            q=q.split(';')[0]+str(i)+':0.5;'
    #print(save)
    if save:
        i=0
        while(i<len(q)):
            if q[i] in ['(',')',','] and q[i+1] not in ['(',')','.']:
                if q.index(':',i)!=-1:
                    j=q.index(':',i)
                    k=q[i+1:j]
                    #print(k)
                    q=q[:i+1]+k.replace(k,label[k],1)+q[j:]
                    i=i+len(k)
            i=i+1
    #print(q)
    return q,label

def input_func(ch,content,gv,gr=False):
    global path_matrix,path_output
    #print("done",gv)
    label={}
    if ch=='c' and gr:       
        dot=graph_to_dot(content)
        content,label=dot_to_newick(dot)
    elif gv:
        #print("done")
        content,label=dot_to_newick(content)
    #print(content)
    handle=io.StringIO(content)
    tree=Phylo.read(handle,"newick")
    Phylo.draw(tree,do_show=False)
    plt.savefig(path_output+'Original_tree.png')
    if path_matrix=="":
        path_matrix=input("Enter the path of cell mutation matrix file.\n")
    matrix=pd.read_csv(path_matrix,header=None,sep=" ")
    matrix=matrix.replace(2,1)
    return tree,matrix,label

def write(content):
    file=open(path_output+'converted_newick.txt','w')
    file.write(content)
    file.close()
    
def phylo_tree(content,typ=0,z=False,gv=False,a=0,call=0,path=None):
    global path_matrix
    if call==1:path_matrix=path
    in_tree,matrix,label=input_func('p',content,gv,call)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 2 for Clonal Tree \n 3 for Mutation Tree \n"))
    if typ==3:
        obj=phylo.Phylo_to_Mut(in_tree,matrix,label,z,gv)
    if typ==2:
        obj=phylo.Phylo_to_Clonal(in_tree,matrix,label,z,gv,a)
    dot=obj.convert()
    #print(dot)
    return dot
    
def clonal_tree(content,typ=0,z=False,gv=False,gr=False,call=0,path=None):
    global path_matrix
    if call==1:path_matrix=path
    in_tree,matrix,label=input_func('c',content,gv,gr,call)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 1 for Phylogenetc Tree \n 3 for Mutation Tree \n"))
    if typ==1:
        obj=clonal.Clonal_to_Phylo(in_tree,matrix,label,gr,z,gv)
    if typ==3:
        obj=clonal.Clonal_to_Mut(in_tree,matrix,label,gr,z,gv)
    dot=obj.convert()
    return dot
    
def muta_tree(content,typ=0,z=False,gv=False,a=0,call=0,path=None):
    global path_matrix
    if call==1:path_matrix=path
    in_tree,matrix,label=input_func('m',content,gv,call)
    if typ==0:
        typ=int(input("Select the tree for conversion:- \n 1 for Phylogenetic Tree \n 2 for Clonal Tree \n"))
    if typ==2:
        obj=mut.Mut_to_Clonal(in_tree,matrix,label,z,gv,a)
    if typ==1:
        obj=mut.Mut_to_Phylo(in_tree,matrix,label,z,gv)
    dot=obj.convert()
    return dot
def main():
    global path_matrix,path_output
    path_input=""
    filename="Converted_Tree"
    typ=0
    dot=""
    gr=False
    z=False
    a=0
    arg=sys.argv
    if len(arg)==1:
        ch=int(input("Enter the tree format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\n"))
        gr=bool(input("The input tree is in igraph format? 1 or 0?\n"))
        z=bool(input("The numbering of cell/mutations in the current tree starts from 0 or 1?\n"))
        a=int(input("Enter the clustering algorithm to be used"))
        
    else:
        arg=arg[1:]
        short_options="t:i:m:o:f:r:gza:"
        long_options=["tree","input","matrix","output","filename","result","igraph","zero","algo"]
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
            if curr_arg in ("-a","--algo"):
                a=int(curr_val)
    #print(gr)
    gv=False
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
    if '.gv' in path_input:
        gv=True 
    content=open(path_input,"r").read()
    #print(content)
    if ch==1:
        dot=phylo_tree(content,typ,z,gv,a)
    elif ch==2:
        dot=clonal_tree(content,typ,z,gv,gr)
    elif ch==3:
        dot=muta_tree(content,typ,z,gv,a)
    else:
        print("Wrong Tree Choice inputted.")
        sys.exit(2)
        
    newick,_=dot_to_newick(dot,save=True)
    dot.render(filename=filename,directory=path_output)
    print("Output will be stored in results folder.")
    write(newick)
if __name__=="__main__":
    main()