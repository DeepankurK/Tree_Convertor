from graphviz import Digraph
import pandas as pd
from collections import Counter 

class Node_2:
    def __init__(self,nodes=[],self_mut=[],br_mut=[],node_label=[]):
        self.nodes=nodes
        self.self_mut=self_mut
        self.node_label=node_label
        
    def show(self):
        ls=[]
        for i in self.nodes:
            ls=ls+i.self_mut
           # print(i.self_mut,"bef")
        if ls==[]: return
        print(self.self_mut,len(self.nodes),ls,"bef")
        ele,cnt=Counter(ls).most_common(1)[0] 
        if cnt>=2:
            new_node=Node_2(br_mut=[],nodes=[],node_label=self.node_label)
            curr_nodes=self.nodes.copy()
            new_node.self_mut.append(ele)
            for j in curr_nodes:
                if ele in j.self_mut:
                    self.nodes.remove(j)
                    j.self_mut.remove(ele)
                    new_node.nodes.append(j)
            print(self.self_mut,len(self.nodes),"now")
            self.nodes.append(new_node)
            self.show()
        else: 
            print(self.self_mut,len(self.nodes),"done")
            for i in self.nodes:
                i.show()    
                
    def optimize(self,root,node):
        node.self_mut=self.node_label[root]
        print(node,node.self_mut)
        for i in root.clades:
            node_curr=Node_2(br_mut=[],nodes=[],node_label=self.node_label)
            node.nodes.append(node_curr)
            if len(i.clades)!=0:
                ret=self.optimize(i,node_curr)
                if ret is not None:
                    for j in ret:
                        node.nodes.append(j)
            else:
                node_curr.self_mut=self.node_label[i]
                print(node_curr,node_curr.self_mut)
        if len(node.self_mut)==1 and len(root.clades)!=0:
            q=node.nodes
            node.nodes=[]
            return q
            


class Clonal_to_Mut:

    def __init__(self,tree,matrix,label):
        self.tree=tree
        self.matrix=matrix
        self.label=label
        self.c=1
        self.root=list(self.tree.find_clades(terminal=False,order='level'))[0]
        self.dot=Digraph(comment="Clonal Tree",format='png')
        self.dot.attr('node')
        self.dot.node('0','N')
        self.node_label={}
        
    def get_mut(self,ls):
        temp_pd=pd.DataFrame()
        #print(ls)
        for i in ls.split(" "):
            #print(i)
            if i!='':
                temp_pd[i]=self.matrix.loc[:,int(i.strip())-1]
        mut_list=[]
        for i in range(len(temp_pd)):
            d=1
            for j in range(len(temp_pd.columns)):
                if temp_pd.iloc[i,j]==0:
                    d=0
            if d==1:
                mut_list.append(str(i+1))
        return mut_list
    
    def transverse(self,node):
        #print(node.name,node.confidence)
        if len(node.clades)==0:
            temp=self.get_mut(self.label[str(node.name)])
            temp.append('C:'+str(node.name))
            #print(temp)
            return temp
        else:
            mut_clus=[]
            muts=[]
            change=[]
            for i in node.clades:
                mut_clus.append(self.transverse(i))
                muts=muts+mut_clus[-1]
            node_list=self.get_mut(self.label[str(node.confidence)])
            node_list.append('C:'+str(node.confidence))
            #print(node.confidence,node_list,mut_clus,"1")
            fr=dict(Counter(muts))
            for x,y in fr.items():
                if y==len(node.clades) and x not in node_list and y!=1:
                    node_list.append(x)
                    change.append(x)
            temp=mut_clus.copy()
            for i in temp:
                mut_clus.remove(i)
                dif=[k for k in node_list if k not in i]
                for j in dif:
                    if 'C' not in j:
                        i.append('-'+j)
                i=[k for k in i if k not in change and k not in node_list]
                mut_clus.append(i)
                #print(i)
            #print(node.confidence,node_list,mut_clus,"2")
            if len(node.clades)!=1:
                comm=mut_clus[0]
                for i in range(1,len(mut_clus)):
                    comm=[j for j in comm if j in mut_clus[i]]
                for i in comm:
                    node_list.append('+'+i)
                for i,z in zip(mut_clus,node.clades):
                    i=[k for k in i if k not in comm]
                    self.node_label[z]=i
                    #print(i)
            else:self.node_label[node.clades[0]]=i
            if self.root==node:self.node_label[node]=node_list
            else:return node_list
            
    def create_node(self,curr,st):
        self.dot.node(str(self.c),st)
        self.dot.edge(str(curr),str(self.c))
        self.c=self.c+1
        
    def tree_to_dot(self,root,prev2):
        j=root.self_mut
        print("prev",j)
        for i in j:
            if '-' in i and '+' not in i:
                self.create_node(prev2,i)
                prev2=self.c-1
        for i in j:
            if '-' not in i and 'C' not in i and '+' not in i:
                self.create_node(prev2,i)
                prev2=self.c-1
    #print(prev,c)
        for i in j:
            if 'C' in i:
                self.dot.node(str(self.c),i)
                self.dot.edge(str(prev2),str(self.c))
                self.c=self.c+1
        for i in j:
            if '+' in i:
                self.create_node(prev2,i.split('+')[1])
                prev2=self.c-1
                
        for j in root.nodes:
            print("imp",j,j.self_mut)
            self.tree_to_dot(j,prev2)
           # elif len(k.nodes)==0:
            #    self.tree_to_dot(k,el_prev)
    def convert(self):  
        self.transverse(self.root)
        #print(self.label)
        #print()
        #print(self.node_label)
        #print()
        root2=Node_2(nodes=[],self_mut=[],br_mut=[],node_label=self.node_label)
        root2.optimize(self.root,root2)
        root2.show()
        self.tree_to_dot(root2,0)
        return self.dot
