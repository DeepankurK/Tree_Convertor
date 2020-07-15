from graphviz import Digraph

#
class Clonal_to_Phylo:
    c=1
    matrix=[]
    dot=Digraph(comment="Clonal Tree",format='png')
    dot.attr('node')
    a=[]
    def __init__(self,in_tree,matri):
        self.a=in_tree[0]
        self.matrix=matri
      
    def transverse(self,node):
        global dot,label
        if label[node.name]!='':
            self.dot.node(label[node.name],label[node.name])
            self.dot.edge(node.name,label[node.name])
        for i in node.descendants:
            self.dot.node(i.name,"")
            self.dot.edge(node.name,i.name)
            self.transverse(i)
            
    def convert(self,in_tree,matri):  
        self.transverse(self.a)
        return self.dot
    
    
