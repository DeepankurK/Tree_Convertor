from Bio import Phylo
import io
import re
import os
from ete3 import Tree
from graphviz import Digraph
import matplotlib.pyplot as plt
import sys
import getopt
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
path_input1=""
path_input2=""
path_output=""
typ=0
filename="Converted_Tree"
dot=""
content=""
class tree:
    def __init__(self,tree,x):
        self.tree=tree
        self.x=x
        self.label={}
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        
    def show(self,node=None):
        if node==None:
            node=self.root
        print(node,node.name,node.confidence,(len(node.clades)))
        for i in node.clades:
            self.show(i)
    def create_label(self,node=None):
        if node==None:
            node=self.root
        if node.confidence != None: self.label[str(node.confidence)]=node
        elif node.name!=None: self.label[str(node.name)]=node
        for i in node.clades:
            self.create_label(i)
            
    def delete_head(self,node):
        if str(node.confidence) in self.label.keys(): del self.label[str(node.confidence)]
        elif str(node.name) in self.label.keys(): del self.label[str(node.name)]
    
    def create_distmat(self):
        for k,x in enumerate(self.label.keys()):
            self.label[x]=[self.label[x],k]
        #print(self.label)
        self.leafs=len(self.label.keys())
        self.distmat=np.zeros((self.leafs, self.leafs))
        ls=list(self.label.keys())
        print(self.tree.get_path(self.label[ls[1]][0]),self.tree.get_path(self.label[ls[2]][0]))
        for i in range(0,self.leafs-1):
            for j in range(i+1,self.leafs):
                temp1=self.tree.get_path(self.label[ls[i]][0])
                temp2=self.tree.get_path(self.label[ls[j]][0])
                self.distmat[i][j]=len([k for k in temp1+temp2 if (k not in temp1 and k in temp2) or (k in temp1 and k not in temp2) ])
        self.sum=np.sum(self.distmat)
        
        
class compare_tree:
    def __init__(self,tree1,tree2,x,y,typ):
        self.obj1=tree(tree1,x)
        self.obj2=tree(tree2,y)
        self.typ=typ
    def pair(self,obj1,obj2):
        obj1.create_distmat()
        obj2.create_distmat()
        print("Tree1 pairwise cell distance:- "+str(obj1.sum))
        print("Tree2 pairwise shortest-path cell distance:- "+str(obj2.sum))
        l=obj1.leafs
        if obj1.leafs!=obj2.leafs:print("Not equal")
        else: print(obj1.leafs)
        #to be done
        
        print("Overall pairwise cell shortest-path distance:- "+ str(abs(obj1.sum-obj2.sum)))
        print("Normalized pairwise cell shortest-path distance:- "+str((2*abs(obj1.sum-obj2.sum))/(l*(l-1))))

    def compare(self):
        self.obj1.show()
        self.obj1.create_label()
        self.obj1.delete_head(self.obj1.root)
        #print(self.obj1.label)
        self.obj2.show()
        self.obj2.create_label()
        self.obj2.delete_head(self.obj2.root)
        #print(self.obj1.label)
        if self.typ==1:self.pair(self.obj1,self.obj2)
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

def input_func(ch,content,num,gr=False):
    label={}
    if ch=='c' and gr:       
        dot=graph_to_dot(content)
        content,label=dot_to_newick(dot)
    handle=io.StringIO(content)
    tree=Phylo.read(handle,"newick")
    Phylo.draw(tree,do_show=False)
    plt.savefig(path_output+'Tree'+str(num)+'.png')
    if ch=='c':
        return tree,label
    return tree

def phylo_tree(content1,content2,typ,x,y):
    tree1=input_func('p',content1,1)
    tree2=input_func('p',content2,2)
    obj=compare_tree(tree1,tree2,x,y,typ)
    obj.compare()

    
def clonal_tree(content1,content2,typ,p,q,x,y):
    tree1,label1=input_func('c',content1,1,p)
    tree2,label2=input_func('c',content2,2,q)
    obj=compare_tree(tree1,tree2,x,y,typ)
    obj.compare()
    
def muta_tree(content1,content2,typ,x,y):
    tree1=input_func('m',content1,1)
    tree2=input_func('m',content2,2)
    obj=compare_tree(tree1,tree2,x,y,typ)
    obj.compare()

    


p=False
q=False
x=False
y=False
arg=sys.argv
if len(arg)==1:
    ch=int(input("Enter the tree format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\nNote:- If you want to just calculate the metrix for a single tree consider 2nd tree as a copy of 1st and give the details.\n"))
    p=bool(input("The Tree1 tree is in igraph format? 1 or 0?\n"))
    q=bool(input("The Tree2 tree is in igraph format? 1 or 0?\n"))
    x=bool(input("The numbering of cell/mutations in the Tree1 starts from 0 or 1?\n"))
    y=bool(input("The numbering of cell/mutations in the Tree1 starts from 0 or 1?\n"))
    
else:
    arg=arg[1:]
    short_options="t:a:b:pqxym:o:"
    long_options=["tree","tree1","tree2","igraph1","igraph2","zero1","zero2","metric","output"]
    try:
        arguments, values = getopt.getopt(arg, short_options, long_options)
    except getopt.error as err:
        print (str(err))
        sys.exit(2)
    for curr_arg, curr_val in arguments:
        if curr_arg in ("-t","--tree"):
            ch=int(curr_val)
        if curr_arg in ("-a","--tree1"):
            path_input1=curr_val
        if curr_arg in ("-b","--tree2"):
            path_input2=curr_val
        if curr_arg in ("-m","--metric"):
            typ=int(curr_val)
        if curr_arg in ("-p","--igraph1"):
            p=True
        if curr_arg in ("-q","--igraph2"):
            q=True
        if curr_arg in ("-x","--zero1"):
            x=True
        if curr_arg in ("-y","--zero2"):
            y=True
        if curr_arg in ("-o","--output"):
            path_output=curr_val


if path_input1=="":
    path_input1=input("Enter the path of 1st input file.\n")
content1=open(path_input1,"r").read()
if path_input2=="":
    path_input2=input("Enter the path of 2nd input file.\n")
content2=open(path_input2,"r").read()
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
if typ==0:typ=int(input("The number of metric to be used.\n"))
if ch==1:
    phylo_tree(content1,content2,typ,x,y)
elif ch==2:
    clonal_tree(content1,content2,typ,p,q,x,y)
elif ch==3:
    muta_tree(content1,content2,typ,x,y)
else:
    print("Wrong Tree Choice inputted.")
    sys.exit(2)
