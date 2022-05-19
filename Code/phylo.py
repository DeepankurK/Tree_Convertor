import pandas as pd
import re
from graphviz import Digraph
from collections import Counter 
from sklearn_extra.cluster import KMedoids
from sklearn import metrics
import numpy as np
import cluster
#back mutation has to be implmented in Phylo_to_mut
#2-21 issue in Phylo_to_clonal


class Node:
    def __init__(self, a='',b='',comm='',n1=None,n2=None):
        self.a=a
        self.b=b
        self.comm=comm
        self.n1=n1
        self.n2=n2  
        
class Node_2:
    def __init__(self,nodes=[],br_mut=[]):
        self.nodes=nodes
        self.br_mut=br_mut
                   
    def show(self):
        ls=[]
        if self.br_mut==[] or None in self.br_mut:return
        for i,j in zip(self.br_mut,self.nodes):
            #print(i,"issue")
            ls=ls+i
        if ls==[] or [] in self.br_mut:print(self.br_mut,"error")
        #print(self.br_mut)
        ele,cnt=Counter(ls).most_common(1)[0] 
        if cnt>=2:
            new_node=Node_2(br_mut=[],nodes=[])
            #print(self.br_mut,len(self.br_mut))
            curr_mut=self.br_mut.copy()
            curr_nodes=self.nodes.copy()
            for i,j in zip(curr_mut,curr_nodes):
                if ele in i:
                    self.br_mut.remove(i)
                    self.nodes.remove(j)
                    i.remove(ele)
                    if i==[]:
                        for k,l in zip(j.br_mut,j.nodes):
                            new_node.br_mut.append(k)
                            new_node.nodes.append(l)
                    else:
                        new_node.br_mut.append(i)
                        new_node.nodes.append(j)

            self.br_mut.append([ele])
            self.nodes.append(new_node)
            #print(self.br_mut,new_node.br_mut)
            #for i in self.br_mut:
             #   print(i)
            self.show()
        else: 
            #print("done2")
            #print(self.br_mut,len(self.nodes))
            for i in self.nodes:
                i.show()    
                
    def optimize(self,root,node):
        a_dif=root.a.split(" ")
        b_dif=root.b.split(" ")
        #comm=root.comm.split(" ")
        #print(comm,a_dif,b_dif)
        #print(node.br_mut)
        if len(a_dif)!=0 and a_dif!=['']:
            node.br_mut.append(a_dif)
            node_curr=Node_2(br_mut=[],nodes=[])
            #print(node_curr.br_mut,"current")
            node.nodes.append(node_curr)
            if root.n1!=None:
                self.optimize(root.n1,node_curr)
        elif root.n1!=None:
            self.optimize(root.n1,node)
        if len(b_dif)!=0 and b_dif!=['']:
            node.br_mut.append(b_dif)
            node_curr=Node_2(br_mut=[],nodes=[])
            #print(node_curr.br_mut,"current")
            node.nodes.append(node_curr)
            if root.n2!=None:
                self.optimize(root.n2,node_curr)
        elif root.n2!=None:
            self.optimize(root.n2,node)
                
class Phylo_to_Mut:
    def __init__(self,in_tree,matrix,label={},z=False,gv=False):
        self.matrix=matrix
        self.tree=in_tree
        self.label=label
        self.gv=gv
        self.c=1
        self.mapp={}
        self.dot=Digraph(comment="Mutation Tree",format='png')
        self.dot.attr('node',shape='circle')
        self.dot.node('0',"N")
        self.z=z
        
    def tree_recursion(self,a,n1,b,n2):
        a_dif = [i for i in a + b if i not in b ]
        b_dif = [i for i in a + b if i not in a ]
        comm=list(set([i for i in a+b if i in a and i in b]))
        a_dif.sort( key = lambda x: (len (x), x))
        b_dif.sort( key = lambda x: (len (x), x))
        comm.sort( key = lambda x: (len (x), x))
        #print(a_dif,b_dif,comm,"arrays")
        #print(n1,n2)
        tr=Node(comm=" ".join(comm),a=" ".join(a_dif),b=" ".join(b_dif),n1=n1,n2=n2)
        return comm,tr
    
    def get_mut(self,node_des):
        temp_pd=pd.DataFrame()
        #print(node_des.name,"desc")
        if self.gv:
            cell=int(re.sub('\D', '',str(self.label[str(node_des.name)])))
        else:
            cell=int(re.sub('\D', '',str(node_des.name)))
        if self.z:temp_pd[cell]=self.matrix.iloc[:,cell]
        else:temp_pd[cell]=self.matrix.iloc[:,cell-1]
        mut_list=[]
        for i in range(len(temp_pd)):
            if temp_pd.loc[i,cell]==1:
                mut_list.append(str(i+1))
        mut_list.append("Cell: "+str(cell))
        #print(mut_list,"mut_list")
        return mut_list
    
    def transverse(self,node):
        #print(node.clades,"node")
        cnt=0
        for j,q in enumerate(node.clades):
            if len(q.clades)!=0:
                cnt=cnt+j+1
        if cnt==0:
            a_=self.get_mut(node.clades[0])
            b_=self.get_mut(node.clades[1])
            return self.tree_recursion(a_,None,b_,None)
        
        elif cnt!=3:
            a_=self.get_mut(node.clades[3-cnt-1])
            b_,n2=self.transverse(node.clades[cnt-1])
            return self.tree_recursion(a_,None,b_,n2)
        else:
            a_,n1=self.transverse(node.clades[0])
            b_,n2=self.transverse(node.clades[1])
            return self.tree_recursion(a_,n1,b_,n2)
    
    def create_node(self,n):
        #print(mapp,n,str(dot).count("[label="+str(n)+"]"))
        if n in self.mapp.keys():
            self.mapp[n][self.c]=self.c
        else: self.mapp[n]={self.c:self.c}
        self.dot.node(str(self.c),str(n))
        self.c=self.c+1
        
    def value_node(self,n,bound,level):
        #print(self.mapp,n,level)
        if str(n)=="0":return "0"
        else:
            for i in range(self.c,self.c-level-1,-1): 
                if i in self.mapp[n].keys():             
                    return self.mapp[n][i]
            for i in range(bound,0,-1):
                if i in self.mapp[n].keys():             
                    return self.mapp[n][i]
    
    def create_edge(self,n,m,bound,level):
        m=self.value_node(m,bound,level)
        n=self.value_node(n,bound,level)
        self.dot.edge(str(n),str(m))
    
    def tree_to_dot(self,root,el_prev):
        #print(root.br_mut,el_prev)
        for i,k in zip(root.br_mut,root.nodes):
            d=self.c
            #print(i,k.br_mut,len(k.nodes))
            if len(i)!=0:
                self.create_node(i[0])
                self.create_edge(el_prev,i[0],d,len(i))
                for j in range(1,len(i)):
                    self.create_node(i[j])
                    self.create_edge(i[j-1],i[j],d,len(i))
                if len(k.nodes)!=0:
                    self.tree_to_dot(k,i[-1])
            elif len(k.nodes)==0:
                self.tree_to_dot(k,el_prev)
    def back_mut(self,root):
        a_dif=root.a.split(" ")
        b_dif=root.b.split(" ")
        if len(a_dif)!=0 and a_dif!=[''] and root.n2!=None:
            a_test=root.n2.a.split(" ")
            b_test=root.n2.b.split(" ")
            for i in a_dif:
                if i in a_test and i not in b_test:
                    pass
                    #print(i,a_dif,a_test,b_test)
                if i in b_test and i not in a_test:
                    pass
                    #print(i,a_dif,a_test,b_test)
        if len(b_dif)!=0 and b_dif!=[''] and root.n1!=None:
            a_test=root.n1.a.split(" ")
            b_test=root.n1.b.split(" ")
            for i in b_dif:
                if i in a_test and i not in b_test:
                    pass
                    #print(i,b_dif,a_test,b_test)
                if i in b_test and i not in a_test:
                    pass
                    #print(i,b_dif,a_test,b_test)
        if root.n1!=None:self.back_mut(root.n1)            
        if root.n2!=None:self.back_mut(root.n2)
        
    def convert(self,typ=False):
        comm,root=self.transverse(self.tree.get_nonterminals()[0])
        #print(root.comm)
        root2=Node_2()
        root2.optimize(root,root2)
        root2.show()
        el_prev=0
        if len(comm)!=0:
            self.create_node(comm[0])
            self.create_edge(0,comm[0],0,len(comm))
            for j in range(1,len(comm)):
                self.create_node(comm[j])
                self.create_edge(comm[j-1],comm[j],0,len(comm))
            el_prev=comm[-1]
       #if typ:self.optimize(root)
        self.tree_to_dot(root2,el_prev)        
        return self.dot

class Phylo_to_Clonal:
    def __init__(self,in_tree,matrix,label={},z=False,gv=False,a=0):
        self.matrix=matrix
        self.tree=in_tree
        self.label={}
        self.gv=gv
        self.dot=Digraph(comment="Clonal Tree",format='png')
        self.dot.attr('node')
        self.dot.node('0','N')
        self.leafs=0
        self.labels=0
        self.cell_labels=[]
        self.allclades=[]
        self.dic_num={}
        self.c=1
        self.cell=min(matrix.shape[1],21)
        self.cell_root={}
        self.z=z
        self.num=0
        self.a=a
        self.root_dic={}
        
    def dist(self):
         self.allclades=list(self.tree.find_clades(terminal=True,order='level'))
         self.leafs=len(self.allclades)             
         distmat=np.zeros((self.leafs, self.leafs))
         for i in range(0,self.leafs-1):
             for j in range(i+1,self.leafs):
                 distmat[i][j]=self.tree.distance(self.allclades[i],self.allclades[j])
         distmat += distmat.transpose()
         return distmat
     
    def create_labels(self,distmat):
        if self.a==1: self.labels=cluster.affinity(distmat)
        elif self.a==2: self.labels=cluster.agglo(distmat,self.cell)
        elif self.a==3: self.labels=cluster.birch(distmat,self.cell)
        elif self.a==4: self.labels=cluster.kmeans(distmat,self.cell)
        elif self.a==5: self.labels=cluster.kmediods(distmat,self.cell)
        elif self.a==6: self.labels=cluster.spectral(distmat,self.cell) 
        #print(self.labels)
        
    def create_cell_labels(self):
        for i in range(0,max(self.labels)+1):
            temp=[]
            for j in range(0,self.leafs):
                if self.labels[j]==i:
                    if self.gv:
                        temp.append(self.label[str(self.allclades[j].name)])
                    else:
                        temp.append(self.allclades[j].name)
            self.cell_labels.append(sorted(temp))
            #print(self.ell_labels)
    def fill_dic(self):
        for i in range(0,max(self.labels)+1):
            self.dic_num[i]=np.count_nonzero(self.labels==i)
                
    def get_root(self):
        self.root_num={}
        for i in range(0,max(self.labels)+1):
            temp=0
            for j in range(0,self.leafs):
                if self.labels[j]==i:
                    if temp==0:
                        temp=self.allclades[j]
                    else:
                        temp=self.tree.common_ancestor(self.allclades[j],temp)
                        #print(self.allclades[j],self.tree.get_path(self.allclades[j]))
                    #print(temp.branch_length,i,self.tree.get_path(temp))
            self.dic_num[i]=temp
            #print(self.dic_num)
            for j  in self.tree.get_path(temp):
                if j not in self.root_num.keys():
                    self.root_num[j]=1
                else:
                    self.root_num[j]+=1
        #print(self.root_num)
        #print(self.dic_num)
        for i in list(self.dic_num.items()):
            #print(i[0],temp)
            a=self.tree.get_path(i[1])
            a.reverse()
            #print(self.oneclade,len(a))
            if a!=[]:
                temp=a[0]
            else:
                temp=self.oneclade
            #print(a,temp.branch_length)
            for j in a:
                #print(self.root_num[j])
                if self.root_num[j]!=1:
                    temp=j
                    break
           # print(temp.branch_length)
            self.dic_num[i[0]]=temp
            #print(i[0],temp)
            #print(self.tree.get_path(temp),i[0])
        #print(":hi",self.dic_num)
        #print(":hi",self.root_num)
    def check(self,ls):
        for j,i in enumerate(ls):
            ls[j]=re.sub('\D', '',str(i))
        return ls
    
    def node(self,add,pr):
        self.dot.node(str(self.c),'C:'+str(self.num))
        self.num=self.num+1
        self.dot.edge(str(pr),str(self.c))
        pr=self.c
        for j in self.check(self.cell_labels[add]):
            self.c=self.c+1
            self.dot.node(str(self.c),'Cell_'+str(j))
            self.dot.edge(str(pr),str(self.c))
        self.c=self.c+1
    
    def assign(self,node,prev,clus):
        #print(self.tree.get_path(node))
        for i in self.tree.get_path(node):
            if self.root_num[i]!=1 and i in self.dic_num.values():
                if i not in self.root_dic.keys() :
                    self.root_dic[i]=[self.c,prev]
                    self.dot.node(str(self.c),' ')
                    self.c=self.c+1
                    self.dot.edge(str(self.root_dic[i][1]),str(self.root_dic[i][0]))
                prev=self.root_dic[i][0]
        if node not in self.root_dic.keys():
            self.node(clus,0)
        else:self.node(clus,self.root_dic[node][0])
            
    def tree_to_dot(self,root,prev_ele):
        #print(root,flag)
        if root==self.oneclade and root in self.dic_num.values():
            for i in self.dic_num.keys():
                if self.dic_num[i]==root:
                    self.node(i,prev_ele)
                #print(self.cell_labels[i])
        if len(root.clades)==0:return
        ls={}
        for i in list(self.dic_num.items()):
            ls[i[0]]=len(self.tree.get_path(i[1]))
        #print(ls)
        ls=sorted(ls.items(), key = lambda kv:(int(kv[1]), int(kv[0])))
        #print(ls)
        #print("don",self.dic_num)
        for i in ls:
            self.assign(self.dic_num[i[0]],0,i[0])
                
    def convert(self):
        distmat=self.dist()
        self.create_labels(distmat)
        self.oneclade=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.create_cell_labels()
        #print(self.labels)
        #print(self.cell_labels,len(self.cell_labels))
        self.fill_dic()
        #print(self.dic_num)
        self.get_root() 
        #print(self.dic_num)
        self.tree_to_dot(self.oneclade,0)
        #print(self.dot)
        return self.dot


