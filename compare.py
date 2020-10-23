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
import convertor
from sklearn.metrics.cluster import *
np.set_printoptions(threshold=sys.maxsize)


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
            table.append((q[0].strip(),q[1].strip(),0.5))
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

def newick_wrapper(content, label):
    i=0
    while(i<len(content)):
        if content[i] in ['(',')',','] and content[i+1] not in ['(',')','.']:
            if content.index(':',i)!=-1:
                j=content.index(':',i)
                k=content[i+1:j]
                #print(k)
                content=content[:i+1]+k.replace(k,label[k],1)+content[j:]
                i=i+len(k)
        i=i+1
    return content

class tree:
    def __init__(self,tree,mapp,z,typ):
        self.tree=tree
        self.z=z
        self.mapp=mapp
        self.labels={}
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.flag=False
        self.typ=typ

    def check(self,node=None):
       # print(node,node.name,node.confidence,(len(node.clades)))
        if self.mapp!={}:
            self.flag=True
    
    def correct(self,st):
        if self.z:return str(int(re.sub('\D','',st))+1)
        else: return re.sub('\D','',st)
        
    def create_label(self,node):
        if self.typ==1 and len(node.clades)==0 and self.mapp!={}:
            if node.confidence !=None: self.labels[self.correct(str(node.confidence))]=node
            elif node.name!=None: self.labels[self.correct(str(node.name))]=node
        elif self.typ==1 and len(node.clades)==0:
            if node.confidence !=None: self.labels[self.correct(self.mapp[str(node.confidence)])]=node
            elif node.name!=None: self.labels[self.correct(self.mapp[str(node.name)])]=node
        if self.typ==2 and self.mapp!={}:
            #print(node.name,node.confidence)
            if (node.name!=None and 'C' in str(self.mapp[str(node.name)]) and 'Ce' not in str(self.mapp[str(node.name)])) or (node.confidence!=None and 'C' in str(self.mapp[str(node.confidence)]) and 'Ce' not in str(self.mapp[str(node.confidence)])):
                for i in node.clades:
                    if i.confidence !=None: self.labels[self.correct(self.mapp[str(i.confidence)])]=node
                    elif i.name!=None: self.labels[self.correct(self.mapp[str(i.name)])]=node
        elif self.typ==2 and node.name!=None and str('C' in node.name) and 'Ce' not in str(node.name):
            for i in node.clades:
                #print(i.name,node,self.correct(str(i.name)))
                if i.confidence !=None and 'C:' not in str(i.confidence): self.labels[self.correct(str(i.confidence))]=node
                elif i.name!=None and 'C:' not in str(i.name): self.labels[self.correct(str(i.name))]=node
        elif self.typ==2 and node.confidence!=None and 'C' in str(node.confidence) and 'Ce' not in str(node.confidence):
            for i in node.clades:
                #print(i.name,node,self.correct(str(i.name)))
                if i.confidence !=None and 'C:' not in str(i.confidence): self.labels[self.correct(str(i.confidence))]=node
                elif i.name!=None and 'C:' not in str(i.name): self.labels[self.correct(str(i.name))]=node
        for i in node.clades:
            self.create_label(i)
            
    def delete_head(self,node):
        if str(node.confidence) in self.labels.keys(): del self.labels[str(node.confidence)]
        elif str(node.name) in self.labels.keys(): del self.labels[str(node.name)]
        if '' in self.labels.keys(): del self.labels['']
        self.set_leafs()
        
    def set_leafs(self):
        self.leafs=len(self.labels.keys())
        
    def create_distmat(self):
        temp=self.labels
        for k,x in enumerate(self.labels.keys()):
            self.labels[x]=[self.labels[x],k]
        self.distmat=np.zeros((self.leafs, self.leafs))
        ls=list(self.labels.keys())
        for i in range(0,self.leafs-1):
            for j in range(i+1,self.leafs):
                temp1=self.tree.get_path(self.labels[ls[i]][0])
                temp2=self.tree.get_path(self.labels[ls[j]][0])
                self.distmat[i][j]=len([k for k in temp1+temp2 if (k not in temp1 and k in temp2) or (k in temp1 and k not in temp2) ])
        self.sum=np.sum(self.distmat)
        self.labels=temp
        
    def create_cluster(self):
        self.clus_label=[0]*self.leafs
        ls=sorted(list(self.labels.keys()))
        self.clus_num={}
        self.num=0
        for i in range(0,self.leafs-1):
            if self.labels[ls[i]] not in self.clus_num.keys():
                self.clus_num[self.labels[ls[i]]]=self.num
                self.num+=1
            self.clus_label[i]=self.clus_num[self.labels[ls[i]]]
        
class compare_tree:
    def __init__(self,tree1='',label1={},tree2='',label2={},m=0,z1=False,z2=False,typ=0):
        self.obj1=tree(tree1,label1,z1,typ)
        self.obj2=tree(tree2,label2,z2,typ)
        self.typ=typ
        self.m=m
        
    def pair(self,obj1,obj2):
        obj1.create_distmat()
        obj2.create_distmat()
        #print(obj1.distmat,obj2.distmat)
        print("Total Leafs are "+str(obj1.leafs))
        print("Tree1 pairwise shortest-path cell distances:- "+str(obj1.sum))
        print("Tree2 pairwise shortest-path cell distances:- "+str(obj2.sum))
        print("Overall pairwise cell shortest-path distances:- "+ str(abs(obj1.sum-obj2.sum)))
        print("Normalized pairwise cell shortest-path distances scorez:- "+str((2*abs(obj1.sum-obj2.sum))/(obj1.leafs*(obj1.leafs-1))))
    
    
    '''
    def rf(self,obj1,obj2):
        Phylo.write(obj1.tree, "temp.txt", "newick")
        content1=open(os.getcwd()+'/temp.txt',"r").read()
        print(content1)
        content1=newick_wrapper(content1,obj1.labels)
        print(content1)
        t1=Tree(content1)
        os.remove(os.getcwd()+'/temp.newick')
        Phylo.write(obj2.tree, "temp.newick", "newick")
        content2=open(os.getcwd()+'/temp.newick',"r").read()
        content2=newick_wrapper(content2,obj2.labels)
        t2=Tree(content2)
        os.remove(os.getcwd()+'/temp.newick')

        parts = t1.robinson_foulds(t2)
        print(parts)
        # to  bw wriiten
        print (t1,t2)
        #print ("RF distance is %s over a total of %s" %(rf, max_rf))
        #print ("Partitions in tree2 that were not found in tree1:", parts_t1 - parts_t2)
        #print ("Partitions in tree1 that were not found in tree2:", parts_t2 - parts_t1)
'''
    def mutual(self,obj1,obj2,m):
        obj1.create_cluster()
        obj2.create_cluster()
        #print(obj1.clus_label,obj2.clus_label)
        print("Total Leafs are "+str(obj1.leafs))
        if m==2:
            print("Adjusted Rand Score :-"+ str(adjusted_rand_score(obj1.clus_label,obj2.clus_label)))
        elif m==3:
            print("Adjusted Mutual Info Score :-"+ str(adjusted_mutual_info_score(obj1.clus_label,obj2.clus_label)))
        elif m==4:
            print("Adjusted Rand Score :-"+ str(adjusted_rand_score(obj1.clus_label,obj2.clus_label)))
        elif m==5:
            print("V Measure Score :-"+ str(v_measure_score(obj1.clus_label,obj2.clus_label)))
        elif m==6:
            print("Calinski Harabasz Score :-"+ str(calinski_harabasz_score([[i] for i in range(0,obj1.leafs)],obj1.clus_label))+" , "+str(calinski_harabasz_score([[i] for i in range(0,self.obj1.leafs)],obj2.clus_label)))
        elif m==7:
             print("Davies Bouldin Score :-"+ str(davies_bouldin_score([[i] for i in range(0,obj1.leafs)],obj1.clus_label))+" , "+str(davies_bouldin_score([[i] for i in range(0,self.obj1.leafs)],obj2.clus_label)))
        elif m==8:
            print("Completeness Score :-"+ str(completeness_score(obj1.clus_label,obj2.clus_label)))
        elif m==9:
            print("Fowlkes Mallows Score :-"+ str(fowlkes_mallows_score(obj1.clus_label,obj2.clus_label)))
        elif m==10:
            print("Homogeneity Score :-"+ str(homogeneity_score(obj1.clus_label,obj2.clus_label)))
        elif m==11:
            print("Normalized Mutual Info Score :-"+ str(normalized_mutual_info_score(obj1.clus_label,obj2.clus_label)))
        elif m==12:
            print("Contingency Matrix :- \n"+ str(contingency_matrix(obj1.clus_label,obj2.clus_label)))
        
    def compare(self):
        self.obj1.check(self.obj1.root)
        self.obj1.create_label(self.obj1.root)
        self.obj1.delete_head(self.obj1.root)
        #print(self.obj1.labels)
        
        self.obj2.check(self.obj1.root)
        self.obj2.create_label(self.obj2.root)
        self.obj2.delete_head(self.obj2.root)
        #print(self.obj2.labels)
        if self.obj1.leafs!=self.obj2.leafs:
            print("Number of leafs are not equal in both trees.\n"+str(self.obj1.leafs)+" in 1st Tree and "+str(self.obj2.leafs)+" in 2nd Tree.")
            temp1=list(self.obj1.labels.keys())
            temp2=list(self.obj2.labels.keys())
            #print(temp1,temp2)
            ls=[k for k in temp1+temp2 if (k not in temp1 and k in temp2) or (k in temp1 and k not in temp2) ]
            for i in ls:
                print("Leaf "+i+" is not present in both trees therefore removed from calculations")
                if i in temp1:del self.obj1.labels[i]
                if i in temp2:del self.obj2.labels[i]
        self.obj1.set_leafs()
        self.obj2.set_leafs()
        if self.m==1:self.pair(self.obj1,self.obj2)
        elif self.m<=12:self.mutual(self.obj1,self.obj2,self.m)
        else: print("Wrong metric choice" )



path_input1=""
path_input2=""
path_output=""
ch1=0
ch2=0
typ=0
dot=""
gr1=False
gr2=False
z1=False
z2=False
gv1=False
gv2=False
r=0
label1={}
label2={}
arg=sys.argv
a=0
path_matrix=""
#print(arg)
if len(arg)==1:
    ch1=int(input("Enter the tree1 format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\nNote:- If you want to just calculate the matrix for a single tree consider 2nd tree as a copy of 1st and give the details.\n"))
    ch2=int(input("Enter the tree1 format for input:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree.\nNote:- If you want to just calculate the matrix for a single tree consider 2nd tree as a copy of 1st and give the details.\n"))
    gr1=bool(input("The Tree1 tree is in igraph format? 1 or 0?\n"))
    gr2=bool(input("The Tree2 tree is in igraph format? 1 or 0?\n"))
    z1=bool(input("The numbering of cell/mutations in the Tree1 starts from 0 or 1?\n"))
    z2=bool(input("The numbering of cell/mutations in the Tree1 starts from 0 or 1?\n"))
    path_matrix=input("Enter the path of cell mutation matrix file.\n")
    r=int(input("Enter the resulting format for output:-\n1. Phylogenetic Tree\n2. Clonal Tree\n3. Mutation Tree."))
else:
    arg=arg[1:]
    short_options="t:u:a:b:pqxym:o:r:l:s:"
    long_options=["choice1","choice2","tree1","tree2","igraph1","igraph2","zero1","zero2","metric","output","result","algo","matrix"]
    try:
        arguments, values = getopt.getopt(arg, short_options, long_options)
    except getopt.error as err:
        print (str(err))
        sys.exit(2)
    for curr_arg, curr_val in arguments:
        if curr_arg in ("-t","--choice1"):
            ch1=int(curr_val)
        if curr_arg in ("-u","--choice2"):
            ch2=int(curr_val)
        if curr_arg in ("-a","--tree1"):
            path_input1=curr_val
        if curr_arg in ("-b","--tree2"):
            path_input2=curr_val
        if curr_arg in ("-m","--metric"):
            typ=int(curr_val)
        if curr_arg in ("-p","--igraph1"):
            gr1=True
        if curr_arg in ("-q","--igraph2"):
            gr2=True
        if curr_arg in ("-x","--zero1"):
            z1=True
        if curr_arg in ("-y","--zero2"):
            z2=True
        if curr_arg in ("-o","--output"):
            path_output=curr_val
        if curr_arg in ("-r","--result"):
            r=int(curr_val)
        if curr_arg in ("-l","--algo"):
            a=int(curr_val)
        if curr_arg in ("-s","--matrix"):
            path_matrix=curr_val


if path_input1=="":
    path_input1=input("Enter the path of 1st input file.\n")
if '.gv' in path_input1:
    gv1=True
content1=open(path_input1,"r").read()
if path_input2=="":
    path_input2=input("Enter the path of 2nd input file.\n")
if '.gv' in path_input2:
    gv2=True
content2=open(path_input2,"r").read()
if path_matrix=="" and (r!=ch1 or r!=ch2):
    path_matrix=input("Enter the path of cell mutation matrix file here.\n")
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

if r==1:
    if ch1==2:
        print("This conversion is not possible")
        sys.exit(2)
    elif ch1==3:
        dot=convertor.muta_tree(content1,1,z1,gv1,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content1,label1=dot_to_newick(dot,save=False)
        z1=False
        gv1=False
    elif gv1==True:
        content1,label1=dot_to_newick(content1)
    if ch2==2:
        print("This conversion is not possible")
        sys.exit(2)
    elif ch2==3:
        dot=convertor.muta_tree(content2,1,z2,gv2,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content2,label2=dot_to_newick(dot,save=False)
        z2=False
        gv2=False
    elif gv2==True:
        content2,label2=dot_to_newick(content2)
elif r==2:
    if ch1==1:
        dot=convertor.phylo_tree(content1,2,z1,gv1,a,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content1,label1=dot_to_newick(dot,save=False)
        z1=False
        gv1=False
    elif ch1==3:
        dot=convertor.muta_tree(content1,2,z1,gv1,a,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content1,label1=dot_to_newick(dot,save=False)
        z1=False
        gv1=False
    elif gv1==True:
        content1,label1=dot_to_newick(content1)
    elif gr1==True:
        dot=graph_to_dot(content1)
        content1,label1=dot_to_newick(dot)
    if ch2==1:
        dot=convertor.phlyo_tree(content2,2,z2,gv2,a,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content2,label2=dot_to_newick(dot,save=False)
        z2=False
        gv2=False
    elif ch2==3:
        dot=convertor.muta_tree(content2,2,z2,gv2,a,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content2,label2=dot_to_newick(dot,save=False)
        z2=False
        gv2=False
    elif gv2==True:
        content2,label2=dot_to_newick(content2)
    elif gr2==True:
        dot=graph_to_dot(content2)
        content2,label2=dot_to_newick(dot)
elif r==3:
    if ch1==1:
        dot=convertor.phylo_tree(content1,3,z1,gv1,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content1,label1=dot_to_newick(dot,save=False)
        z1=False
        gv1=False
    elif ch1==2:
        dot=convertor.clonal_tree(content1,3,z1,gv1,gr1,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content1,label1=dot_to_newick(dot,save=False)
        z1=False
        gv1=False
    elif gv1==True:
        content1,label1=dot_to_newick(content1)
    if ch2==1:
        dot=convertor.phylo_tree(content2,3,z2,gv2,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content2,label2=dot_to_newick(dot,save=False)
        z2=False
        gv2=False
    elif ch2==2:
        dot=convertor.clonal_tree(content2,3,z2,gv2,gr2,1,path_matrix)
        dot.render(filename="Converted_Tree",directory=path_output)
        content2,label2=dot_to_newick(dot,save=False)
        z2=False
        gv2=False
    elif gv2==True:
        content2,label2=dot_to_newick(content2)
else:
    print("Wrong Tree Choice inputted.")
    sys.exit(2)
    
tree1=Phylo.read(io.StringIO(content1),"newick")
Phylo.draw(tree1,do_show=False)
plt.savefig(path_output+'Original_tree1.png')

tree2=Phylo.read(io.StringIO(content2),"newick")
Phylo.draw(tree2,do_show=False)
plt.savefig(path_output+'Original_tree2.png')
#print(content2,label2)
obj=compare_tree(tree1,label1,tree2,label2,typ,z1,z2,r)
ret=obj.compare()