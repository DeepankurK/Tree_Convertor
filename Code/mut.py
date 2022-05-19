import pandas as pd
from graphviz import Digraph
import warnings
warnings.filterwarnings("ignore")
from collections import Counter 
import numpy as np
from sklearn_extra.cluster import KMedoids
from sklearn import metrics
import re
import cluster
#2 to 21 issue to be seen

class Mut_to_Clonal:
    def __init__(self,in_tree,matrix,label={},z=False,gv=False,a=0):
        self.tree=in_tree
        self.matrix=matrix
        self.mapp={}
        self.label=label
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.c=1
        self.dot=Digraph(comment="Clonal Tree",format='png')
        self.dot.attr('node')
        self.dot.node('0','N')
        self.mut_map={}
        self.cell_map={}
        self.labels=[]
        self.dic_num={}
        self.place={}
        self.cell=min(matrix.shape[1],21)
        self.z=z
        self.mut_ls=[]
        self.gv=gv
        self.num=0
        self.a=a
    def mut_mean(self,node_des,prev_mut):
        temp_pd=pd.DataFrame()
        #print(node_des)
        for i in node_des:
            #print(i,i.name,i.confidence)
            if i.confidence!=None:muta=int(i.confidence)
            elif i.name!=None:muta=int(i.name)
            else:muta=len(self.matrix)+1
            if self.label!={}:
                 muta=int(re.sub('\D', '',str(self.label[str(muta)])))
            if int(muta)<=len(self.matrix):
                if self.z:temp_pd[muta]=self.matrix.loc[muta,:]
                else:temp_pd[muta]=self.matrix.loc[muta-1,:]
        mut_list=[]
        for i in range(len(temp_pd)):
            d=1
            for j in range(len(temp_pd.columns)):
                if temp_pd.iloc[i,j]==0:
                    d=0
            if d==1:
                mut_list.append(str(i+1))
        for i in prev_mut:
            mut_list=list(set(mut_list)-set(i))
            mut_list.sort( key = lambda x: (len (x), x))
        #print(mut_list,"mit_list")
        return mut_list
    
        
    def transverse(self,nodes):
        #print(nodes[-1],nodes[-1].confidence,nodes[-1].clades)
        if len(nodes[-1].clades)==0:
            ls=self.mut_mean(nodes,[[]])
            self.cell_map[nodes[-1]]=ls
            for i in ls:
                self.mut_map[i]=nodes[-1]
            return ls
        else:
            mut_clus=[]
            for i in nodes[-1].clades:
                temp_nodes=nodes.copy()
                temp_nodes.append(i)
                mut_clus.append(self.transverse(temp_nodes))
            #print(mut_clus,"done")
            node_list=self.mut_mean(nodes,mut_clus)
            self.cell_map[nodes[-1]]=node_list
            for i in node_list:
                self.mut_map[i]=nodes[-1]
            #print(nodes[-1].confidence,node_list,mut_clus,"help")
            for i in mut_clus:
                node_list=list(set(i+node_list))
            return node_list
                 
    def create_dist(self):
        self.keys=list(self.mut_map.keys())
        self.keys.sort( key = lambda x: (len (x), x))
        self.distmat=np.zeros((len(self.keys),len(self.keys)))
        for i in range(0,len(self.keys)):
            for j in range(0,len(self.keys)):
                    a=self.tree.get_path(self.mut_map[self.keys[i]])
                    b=self.tree.get_path(self.mut_map[self.keys[j]])
                    self.distmat[i][j]=len(list(set([i for i in a+b if (i in a and i not in b) or (i in b and i not in a)])))
                    
    
    def create_labels(self,distmat):
        #print(distmat)
        if self.a==1: self.labels=cluster.affinity(distmat)
        elif self.a==2: self.labels=cluster.agglo(distmat,self.cell)
        elif self.a==3: self.labels=cluster.birch(distmat,self.cell)
        elif self.a==4: self.labels=cluster.kmeans(distmat,self.cell)
        elif self.a==5: self.labels=cluster.kmediods(distmat,self.cell)
        elif self.a==6: self.labels=cluster.spectral(distmat,self.cell)   
        for i in range(len(self.keys)):
            self.mut_map[self.keys[i]]=[self.mut_map[self.keys[i]],self.labels[i]]
    
    def fill_dic(self):
        self.dic_num={}
        for i in range(0,max(self.labels)+1):
            if np.count_nonzero(self.labels==i)!=0:
                self.dic_num[i]=np.count_nonzero(self.labels==i)
                
    
    def set_place(self,root,a_i):
        if a_i not in self.place.keys(): 
            self.place[a_i]=root
        else:
            for i in self.tree.get_path(self.place[a_i]):
                if i in self.tree.get_path(root):
                    self.place[a_i]=root
    def placing(self,root):
        #print(root.confidence,root.name)
        if len(root.clades)==1 and len(self.cell_map[root])==0:
            return self.placing(root.clades[0])
        else:
            for i in root.clades:
               a,a_i=self.placing(i) 
               #print(a,a_i)
               if a_i in self.dic_num.keys():
                   self.dic_num[a_i]-=a
                   self.set_place(root,a_i)
            #print(self.dic_num,root.confidence,root.name)
            if len(self.cell_map[root])>0:
                temp=[]
                for i in self.cell_map[root]:
                    temp.append(self.mut_map[i][1])
                ele,cnt=Counter(temp).most_common(1)[0]
                temp=dict(Counter(temp))
                for i in temp.keys():
                    if i!=ele:
                        self.dic_num[i]-=temp[i]
                        self.set_place(root,i)
                return  cnt,ele
            else: 
                return 0,max(self.labels)+1
    def cluster(self):
        self.labels={}
        for i in self.place.keys():
            temp=[]
            for j in self.mut_map.keys():
                if i==self.mut_map[j][1]:
                    temp.append(j)
            self.labels[i]=temp
    def mut_show(self):
        for i in self.place.keys():
            #print(self.tree.get_path(self.place[i]),i)
            self.mut_ls=list(set(self.mut_ls+self.tree.get_path(self.place[i])))
        #self.mut_ls=[i.confidence if i.confidence!=None else i.name for i in self.mut_ls]
        
    def show(self,root,prev):
        #print(prev,root.confidence,root.clades,"prev")
        if root in self.mut_ls:
            if root.confidence!=None:
                self.dot.node(str(self.c),str(root.confidence))
            else:
                self.dot.node(str(self.c),str(root.name))
            self.dot.edge(str(prev),str(self.c))
            self.c=self.c+1
            prev=self.c-1
        for i in self.place.keys():
            if root==self.place[i]:
                self.dot.node(str(self.c),"C:"+str(self.num))
                self.num=self.num+1
                self.dot.edge(str(prev),str(self.c))
                pr=self.c
                for j in self.labels[i]:
                    self.c=self.c+1
                    self.dot.node(str(self.c),"Cell_"+str(j))
                    self.dot.edge(str(pr),str(self.c))
                self.c=self.c+1
        #print(prev,root.confidence)
        for i in root.clades:
            self.show(i,prev)              
        
    def convert(self):
        self.transverse([self.root])
        #print(self.mut_map)
        #print(self.cell_map)
        self.create_dist()
        self.create_labels(self.distmat)
        #print(self.mut_map)
        self.fill_dic()
        self.placing(self.root)
        #print(self.place)
        self.mut_show()
        #print(self.mut_ls)
        self.cluster()
        #print(self.labels)
        self.show(self.root,0)
        return self.dot
    

class Mut_to_Phylo:
    def __init__(self,in_tree,matrix,label={},z=False,gv=False):
        self.tree=in_tree
        self.matrix=matrix
        self.mapp={}
        self.label=label
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.c=1
        self.dot=Digraph(comment="Clonal Tree",format='png')
        self.dot.attr('node')
        self.dot.node('0','N')
        self.mut_map={}
        self.labels=[]
        self.z=z
        self.gv=gv
   
    def mut_mean(self,node_des,prev_mut):
        temp_pd=pd.DataFrame()
        for i in node_des:
            if i.confidence!=None:muta=int(i.confidence)
            else :muta=int(i.name)
            if self.gv:
                muta=int(re.sub('\D', '',str(self.label[str(muta)])))
            if int(muta)<=len(self.matrix):
                if self.z:temp_pd[muta]=self.matrix.loc[muta,:]
                else: temp_pd[muta]=self.matrix.loc[muta-1,:]
        mut_list=[]
        for i in range(len(temp_pd)):
            d=1
            for j in range(len(temp_pd.columns)):
                if temp_pd.iloc[i,j]==0:
                    d=0
            if d==1:
                mut_list.append(str(i+1))
        for i in prev_mut:
            mut_list=list(set(mut_list)-set(i))
            mut_list.sort( key = lambda x: (len (x), x))
        #print(mut_list,"mit_list")
        return mut_list
    
        
    def transverse(self,nodes):
        #print(nodes[-1],nodes[-1].confidence,nodes[-1].clades)
        if len(nodes[-1].clades)==0:
            ls=self.mut_mean(nodes,[[]])
            self.mut_map[nodes[-1]]=ls
                #print(i,nodes[-1])
            return ls
        else:
            mut_clus=[]
            for i in nodes[-1].clades:
                temp_nodes=nodes.copy()
                temp_nodes.append(i)
                mut_clus.append(self.transverse(temp_nodes))
            #print(mut_clus,"done")
            node_list=self.mut_mean(nodes,mut_clus)
            self.mut_map[nodes[-1]]=node_list
            #print(nodes[-1].confidence,node_list,mut_clus,"help")
            for i in mut_clus:
                node_list=list(set(i+node_list))
            return node_list
    
    def make(self,curr,nex,st):
        self.dot.node(str(nex),str(st))
        self.dot.edge(str(curr),str(nex))
        self.c=self.c+1
        
    def join(self,ls,prev):
        l=len(ls)
        if l==0: return 
        while(l!=1):
            self.make(prev,self.c," ")
            prev=self.c-1
            self.make(prev,self.c,'Cell: '+ls[l-1])
            l-=1
        self.make(prev,self.c,'Cell: '+ls[0])
    def name(self,i):
        st=''
        if self.gv:
            if i.confidence==None:st=str(self.label[str(i.name)])
            else: st=str(self.label[str(i.confidence)])
        else:
            if i.confidence==None:st=str(i.name)
            else: st=str(i.confidence)
        return st
    def reconstruction(self,root,prev):
        if len(root.clades)>1:
            if len(self.mut_map[root])>=1:
                self.join(self.mut_map[root],prev)
                self.make(prev,self.c," ")
                prev=self.c-1
            l=len(root.clades)
            while(l!=2):
                self.make(prev,self.c,self.name(root.clades[l-1]))
                self.reconstruction(root.clades[l-1],self.c-1)  
                self.make(prev,self.c," ")
                prev=self.c-1
                l-=1
            self.make(prev,self.c,self.name(root.clades[1]))
            self.reconstruction(root.clades[1],self.c-1)
            self.make(prev,self.c,self.name(root.clades[0]))
            self.reconstruction(root.clades[0],self.c-1)
        elif  len(root.clades)==1 and len(self.mut_map[root])>=1:
            self.join(self.mut_map[root],prev)
            self.make(prev,self.c,self.name(root.clades[0]))
            self.reconstruction(root.clades[0],self.c-1)
        elif  len(root.clades)==1:
            self.make(prev,self.c,self.name(root.clades[0]))
            self.reconstruction(root.clades[0],self.c-1)
        elif  len(self.mut_map[root])>=1:
            self.join(self.mut_map[root],prev)
            
    def convert(self):
        self.transverse([self.root])
        #print(self.mut_map)
        self.reconstruction(self.root,0)
        return self.dot
    
    
