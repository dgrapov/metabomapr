#manual calculation
#this requires metabomapr R package: 
devtools::install_github('dgrapov/metabomapr')

#setup
library(metabomapr)
library(dplyr)

#input
#expects that workind directory contains the .csv
data<-read.csv('example.csv',header = TRUE)

#chemical similarity
type<-'CID'
id<-data$PUBCHEM %>% unique()
cid_el<-CID_tanimoto(id,as='edge.list')

#save to .csv
write.csv(cid_el,file='manual_edges.csv',row.names = FALSE)


id<-c('C00123','C01685')
type<-'KEGG'
kegg_el<-get_KEGG_edgeList(id)
