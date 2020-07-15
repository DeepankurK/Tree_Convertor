import pandas as pd
import re
from graphviz import Digraph


class Mut_to_Clonal:
    mapp={}
    c=1
    matrix=pd.DataFrame()
    a=0
    dot=Digraph(comment="Clonal Tree",format='png')
    dot.attr('node')
    def __init__(self,in_tree,matri):
        self.a=in_tree[0]
        self.matrix=matri

    def tree_recursion(self,node_clus,child_clus):
        t=[k for k in node_clus if 'c' in str(k)]
        temp_clus=node_clus.copy()
        self.dot.node(t[0],(" , ").join(node_clus))
        for j,i in enumerate(child_clus):
            child_clus[j]=list(set(child_clus[j])-set(node_clus))
            b=[k for k in child_clus[j] if 'c' in str(k)]
            b_not=[k for k in child_clus[j] if 'c' not in str(k)]
            if b[0]+" [label="+b[0]+']' not in str(self.dot):
                self.dot.node(b[0],(" , ").join(i))
            self.dot.edge(t[0],b[0])
            temp_clus=temp_clus+b_not
        return temp_clus
    
    def mut_mean(self,node_des,prev_mut):
        temp_pd=pd.DataFrame()
        for i in node_des:
            if 'c' not in i.name and 'C' not in i.name:
                muta=int(re.sub('\D', '',str(i.name)))
                if int(muta)<=len(self.matrix):
                    muta=int(muta)
                    temp_pd[muta]=self.matrix.loc[muta-1,:]
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
        if len(mut_list)!=0:
            mut_list.append("c"+str(self.c))
            self.c=self.c+1
       # print(mut_list)
        return mut_list
    
        
    def transverse(self,nodes):
        #print(nodes)
        if len(nodes[-1].descendants)==0:
            return self.mut_mean(nodes,[[]])
        elif len(nodes[-1].descendants)==1:
            temp_nodes=nodes.copy()
            temp_nodes.append(nodes[-1].descendants[0])
            return self.transverse(temp_nodes)
        else:
            mut_clus=[]
            for i in nodes[-1].descendants:
                temp_nodes=nodes.copy()
                temp_nodes.append(i)
                mut_clus.append(self.transverse(temp_nodes))
                #print(mut_clus,"done")
            node_list=self.mut_mean(nodes,mut_clus)
            #print(node_list,mut_clus,"help")
            return self.tree_recursion(node_list,mut_clus)
    
    def convert(self):
        self.transverse([self.a])
        return self.dot
    
    
    
