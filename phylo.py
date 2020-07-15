import pandas as pd
import re
from graphviz import Digraph

class Phylo_to_Mut:
    mapp={}
    c=1
    matrix=pd.DataFrame()
    prob=0.90
    dot=Digraph(comment="Mutation Tree",format='png')
    dot.attr('node',shape='circle')
    dot.node('C','C')
    a=0
    def __init__(self,in_tree,matri):
        self.matrix=matri
        self.a=in_tree[0]
    
    def create_node(self,n):
        #print(mapp,n,str(dot).count("[label="+str(n)+"]"))
        if str(self.dot).count("[label="+str(n)+"]")>=1:
            self.mapp[n]=str(n)+"_"+str(str(self.dot).count("[label="+str(n)+"]"))
            self.dot.node(self.mapp[n],str(n))
        elif 'c' not in str(n):
            self.mapp[n]=str(n)
            self.dot.node(str(n),str(n))
        else:
            self.dot.node(str(n),str(n))
    
    def value_node(self,n,map_temp):
        if n in map_temp.keys() and 'c' not in str(n):
            #print(mapp[n])
            return map_temp[n]
        else: return n
    
    def create_edge(self,n,m,map_temp):
        n=self.value_node(n,map_temp)
        m=self.value_node(m,map_temp)
        self.dot.edge(str(n),str(m))
    
    def tree_recursion(self,a_,a_map,b_,b_map):
        a=a_[0]
        b=b_[0]
        a_dif = [i for i in a + b if i not in b ]
        b_dif = [i for i in a + b if i not in a ]
        comm=list(set([i for i in a+b if i in a and i in b]))
        #print(a_dif,b_dif,comm)
        for i in a_dif+b_dif:
            self.create_node(i)
        ret=[]
        ret.append(comm)
        if len(comm)==0:
            if len(a_dif)!=0:
                self.create_edge(0,a_dif[0],a_map)
                for j in range(1,len(a_dif)):
                    self.create_edge(a_dif[j-1],a_dif[j],a_map)
                for i in range(1,len(a_)):
                    self.create_edge(a_dif[-1],a_[i],a_map)
                #ret.append(value_node(a_dif[0],a_map))
            if len(b_dif)!=0:
                self.create_edge(0,b_dif[0],b_map)
                for j in range(1,len(b_dif)):
                    self.create_edge(b_dif[j-1],b_dif[j],b_map)
                for i in range(1,len(b_)):
                    self.create_edge(b_dif[-1],b_[i],b_map)
                #ret.append(value_node(b_dif[0],b_map))
            return ret,self.mapp
        if  len(a_dif)!=0 and str(a_dif[-1])==a_dif[-1]:
            for j in range(1,len(a_dif)):
                self.create_edge(a_dif[j-1],a_dif[j],a_map)
            ret.append(self.value_node(a_dif[0],a_map))
            
        elif len(a_dif)!=0 and len(comm)!=0:
            for i in range(1,len(a_)):
                self.create_edge(a_dif[-1],a_[i],a_map)
            for j in range(1,len(a_dif)):
                self.create_edge(a_dif[j-1],a_dif[j],a_map)
            ret.append(self.value_node(a_dif[0],a_map))
            
        elif len(a_dif)==0:
                for i in range(1,len(a_)):
                    ret.append(a_[i])
                    
        if len(b_dif)!=0 and str(b_dif[-1])==b_dif[-1] :
            for j in range(1,len(b_dif)):
                self.create_edge(b_dif[j-1],b_dif[j],b_map)
            ret.append(self.value_node(b_dif[0],b_map))
            
        elif len(b_dif)!=0 and len(comm)!=0:
            for i in range(1,len(b_)):
                self.create_edge(b_dif[-1],b_[i],b_map)
            for j in range(1,len(b_dif)):
                self.create_edge(b_dif[j-1],b_dif[j],b_map)
            ret.append(self.value_node(b_dif[0],b_map))
            
        elif len(b_dif)==0:
                for i in range(1,len(b_)):
                    ret.append(b_[i])   
                    
        return ret,self.mapp
    
    def mut_mean(self,node_des):
        temp_pd=pd.DataFrame()
        #print(node.descendants)
        for i in node_des:
            muta=int(re.sub('\D', '',str(i.name)))
            #print(muta)
            temp_pd[muta]=self.matrix.loc[:,muta+1]
        mut_list=[]
        temp_pd=temp_pd.mean(axis=1)
        for j in range(len(temp_pd)):
            if temp_pd[j]>=self.prob:
                mut_list.append(j+1)
        mut_list.append("c"+str(self.c))
        self.mapp["c"+str(self.c)]="c"+str(self.c)
        self.c=self.c+1
        return [mut_list],self.mapp
    
    def transverse(self,node):
        #print(node)
        node_0=[]
        for q in node.descendants:
            if len(q.descendants)==0:
                node_0.append(q)
        if len(node_0)==len(node.descendants):
            return self.mut_mean(node.descendants)
        elif len(node_0)!=0:
            a_,a_map=self.mut_mean(node_0)
            b_,b_map=self.transverse(list(set(node.descendants)-set(node_0))[0])
            return self.tree_recursion(a_,a_map,b_,b_map)
        else:
            a_,a_map=self.transverse(node.descendants[0])
            b_,b_map=self.transverse(node.descendants[1])
            return self.tree_recursion(a_,a_map,b_,b_map)
        
    def convert(self):
        self.transverse(self.a)
        return self.dot
    






