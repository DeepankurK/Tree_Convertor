set.seed(1)
library(oncoNEM)
data_temp <- read.table(file.choose(), header = FALSE, sep = " ", dec = ".")
data_temp <- as.matrix(data_temp)
print(data_temp)
#if index is there then just remove it
data_req <-data_temp[,-1] 
dim(data_req)
oNEM <- oncoNEM$new(Data = data_req,
                    FPR = 0.05,
                    FNR = 0.2)
oNEM$search(delta = 50)
oNEM$best
oNEM$best$tree
plotTree(tree = oNEM$best$tree,clones = NULL,vertex.size = 4)
oNEM.expanded <- expandOncoNEM(oNEM,epsilon = 10,delta = 200,
                               checkMax = 10000,app = TRUE)
plotTree(tree = oNEM.expanded$best$tree,clones = NULL,vertex.size = 25)
oncoTree <- clusterOncoNEM(oNEM = oNEM.expanded,
                           epsilon = 10)
post <- oncoNEMposteriors(tree = oncoTree$g,
                          clones = oncoTree$clones,
                          Data = oNEM$Data,
                          FPR = oNEM$FPR,
                          FNR = oNEM$FNR)
edgeLengths = colSums(post$p_theta)[-1]
plotTree(tree = oncoTree$g,
         clones = oncoTree$clones,
         e.length = edgeLengths,
         label.length = 4,
         axis = TRUE)
print(oncoTree)
capture.output(oncoTree, file = "MyNewFile.txt")