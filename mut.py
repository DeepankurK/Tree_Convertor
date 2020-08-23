import pandas as pd
from graphviz import Digraph
import warnings
warnings.filterwarnings("ignore")
from collections import Counter 
import numpy as np
from sklearn_extra.cluster import KMedoids
from sklearn import metrics

#2 to 21 issue to be seen

class Mut_to_Clonal:
    def __init__(self,in_tree,matrix,z):
        self.tree=in_tree
        self.matrix=matrix
        self.mapp={}
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
   
    def mut_mean(self,node_des,prev_mut):
        temp_pd=pd.DataFrame()
        for i in node_des:
            if i.confidence!=None:muta=int(i.confidence)
            else :muta=int(i.name)
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
                    
    def sil_score(self):
        max_score=0
        for i in range(2,self.cell):
            avg_score=0
            kmedoids = KMedoids(n_clusters=i, random_state=0).fit(self.distmat)
            #print(i,kmedoids.labels_)
            if max(kmedoids.labels_)>0:avg_score=metrics.silhouette_score(self.distmat, kmedoids.labels_)
            #print(avg_score)
            if max_score<avg_score:
                max_score=avg_score
                self.labels=kmedoids.labels_
        #print(self.labels)
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
            
    def show(self,root,prev):
        #print(prev,root.confidence,root.clades,"prev")
        for i in self.place.keys():
            if root==self.place[i]:
                self.dot.node(str(self.c)," ".join(self.labels[i]))
                self.dot.edge(str(prev),str(self.c))
                self.c=self.c+1
                prev=self.c-1
        #print(prev,root.confidence)
        for i in root.clades:
            self.show(i,prev)              
        
    def convert(self):
        self.transverse([self.root])
        self.create_dist()
        self.sil_score()
        self.fill_dic()
        self.placing(self.root)
        self.cluster()
        #print(self.place)
        self.show(self.root,0)
        return self.dot
    

class Mut_to_Phylo:
    def __init__(self,in_tree,matrix,z):
        self.tree=in_tree
        self.matrix=matrix
        self.mapp={}
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.c=1
        self.dot=Digraph(comment="Clonal Tree",format='png')
        self.dot.attr('node')
        self.dot.node('0','N')
        self.mut_map={}
        self.labels=[]
        self.z=z

   
    def mut_mean(self,node_des,prev_mut):
        temp_pd=pd.DataFrame()
        for i in node_des:
            if i.confidence!=None:muta=int(i.confidence)
            else :muta=int(i.name)
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
            self.make(prev,self.c,'Cell:'+ls[l-1])
            l-=1
        self.make(prev,self.c,'Cell:'+ls[0])
    def name(self,i):
        if i.confidence==None: return 'Mut:'+str(i.name)
        return 'Mut:'+str(i.confidence)
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
    
    
