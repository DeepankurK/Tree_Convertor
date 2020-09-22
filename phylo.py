import pandas as pd
import re
from graphviz import Digraph
from collections import Counter 
from sklearn_extra.cluster import KMedoids
from sklearn import metrics
import numpy as np

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
        if ls==[] or [] in self.br_mut:print(self.br_mut,"again")
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
    def __init__(self,in_tree,matrix,label={},z=False,gv=False):
        self.matrix=matrix
        self.tree=in_tree
        self.label={}
        self.gv=gv
        self.dot=Digraph(comment="Mutation Tree",format='png')
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
        
    def dist(self):
         self.allclades=list(self.tree.find_clades(terminal=True,order='level'))
         self.leafs=len(self.allclades)             
         distmat=np.zeros((self.leafs, self.leafs))
         for i in range(0,self.leafs-1):
             for j in range(i+1,self.leafs):
                 distmat[i][j]=self.tree.distance(self.allclades[i],self.allclades[j])
         distmat += distmat.transpose()
         return distmat
     
    def sil_score(self,distmat):
        max_score=0
        for i in range(2,self.cell):
            #print(i)
            kmedoids = KMedoids(n_clusters=i, random_state=0).fit(distmat)
            #print(kmedoids)
            avg_score = metrics.silhouette_score(distmat, kmedoids.labels_)
            if max_score<avg_score:
                max_score=avg_score
                self.labels=kmedoids.labels_
    def create_cell_labels(self):
        for i in range(0,max(self.labels)+1):
            temp=[]
            for j in range(0,self.leafs):
                if self.labels[j]==i:
                    if self.gv:
                        temp.append(self.label[str(self.allclades[j].name)])
                    else:
                        temp.append(self.allclades[j].name)
            self.cell_labels.append(temp)
    def fill_dic(self):
        for i in range(0,max(self.labels)+1):
            self.dic_num[i]=np.count_nonzero(self.labels==i)
            
    def show2(self,root):
        if len(root.clades)==0:
            return -1
        for i in root.clades:
            ret=self.show2(i)
            if ret==-1:
                self.cell_root[i.name]=root
                
    def get_root(self):
        for p,i in enumerate(self.cell_labels):
            comm=self.tree.get_path(self.cell_root[i[0]])
            for j in i[1:]:
                ls=self.tree.get_path(self.cell_root[j])
                comm=[k for k in ls if k in comm]
            if len(comm)!=0:
                self.dic_num[p]=comm[-1]
            else : self.dic_num[p]=self.oneclade
        
        
    def optimize(self):
        for i in range(0,max(self.labels)+1):
            if i in self.dic_num.keys():
                if self.dic_num[i] not in self.dic_num.keys():
                    self.dic_num[self.dic_num[i]]=[i] 
                else:
                    self.dic_num[self.dic_num[i]].append(i) 
                #print(self.dic_num[i].name,i)
                del self.dic_num[i]

        self.comm=[]
        #print(self.dic_num,"done")
        for j in self.dic_num.keys():
            self.comm=list(set(self.comm+self.tree.get_path(j)))
            #print(self.tree.get_path(j))
        #print(self.comm)
    def check(self,ls):
        for j,i in enumerate(ls):
            ls[j]=re.sub('\D', '',str(i))
        return ls
    
    def node(self,typ,add,pr):
        if typ==1:
            self.dot.node(str(self.c),'C: '+" ".join(self.check(self.cell_labels[add])))
            self.dot.edge(str(pr),str(self.c))
            self.c=self.c+1
        else:
            self.dot.node(str(self.c)," ")
            self.dot.edge(str(pr),str(self.c))
            self.c=self.c+1
            
    def tree_to_dot(self,root,prev_ele):
        #print(root,flag)
        if root==self.oneclade and root in self.dic_num.keys():
            for i in self.dic_num[root]:
                self.node(1,i,prev_ele)
                #print(self.cell_labels[i])
        if len(root.clades)==0:return
        a,b=root.clades
        #print(a.branch_length,b.branch_length)
        cnt=0
        if a in self.dic_num.keys():cnt+=1
        if b in self.dic_num.keys():cnt+=2
        if cnt==0:
            if a in self.comm: 
                 self.node(0,a,prev_ele)
                 self.tree_to_dot(a,self.c-1)
            if b in self.comm:
                 self.node(0,b,prev_ele)
                 self.tree_to_dot(b,self.c-1)
        elif cnt==1:
            for i in self.dic_num[a]:
                self.node(1,i,prev_ele)
            self.tree_to_dot(a,self.c-1)
            if b in self.comm:
                self.node(0,b,prev_ele)
                self.tree_to_dot(b,self.c-1)
        elif cnt==2:
            for i in self.dic_num[b]:
                self.node(1,i,prev_ele)
            self.tree_to_dot(b,self.c-1)
            if a in self.comm:
                self.node(0,a,prev_ele)
                self.tree_to_dot(a,self.c-1)
        elif cnt==3:
            for i in self.dic_num[a]:
                self.node(1,i,prev_ele)
            self.tree_to_dot(a,self.c-1)
            for i in self.dic_num[b]:
                self.node(1,i,prev_ele)
                #print(i)
            self.tree_to_dot(b,self.c-1)
                
    def convert(self):
        distmat=self.dist()
        self.sil_score(distmat)
        self.oneclade=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.create_cell_labels()
        #print(self.labels)
        #print(self.cell_labels,len(self.cell_labels))
        self.fill_dic()
        #print(self.dic_num)
        self.show2(self.oneclade)
        self.get_root() 
        #print(self.dic_num)
        self.optimize()
        #print(self.dic_num)
        self.tree_to_dot(self.oneclade,0)
        return self.dot


