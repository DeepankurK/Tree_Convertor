from graphviz import Digraph

#
class Clonal_to_Phylo:
    c=1
    matrix=[]
    dot=Digraph(comment="Clonal Tree",format='png')
    dot.attr('node')
    a=[]
    label={}
    def __init__(self,in_tree,matri,label):
        self.a=in_tree[0]
        self.matrix=matri
        self.label=label
      
    def transverse(self,node):
        if self.label[node.name]!='':
            self.dot.node(self.label[node.name],self.label[node.name])
            self.dot.edge(node.name,self.label[node.name])
        for i in node.descendants:
            self.dot.node(i.name,"")
            self.dot.edge(node.name,i.name)
            self.transverse(i)
            
    def convert(self):  
        self.transverse(self.a)
        return self.dot
    
    
