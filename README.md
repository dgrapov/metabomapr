
### Install

```r
#devtools::install_github("dgrapov/metabomapr")
library("metabomapr")
```

### Examples
#### Convert PubChem CID to tanimoto similarity adjacency matrix

```r
cids<-c("[]","51", "440649","[]", "439161",  NA )
(res<-CID_tanimoto(cids))
```

```
##               51    440649   439161
## 51     1.0000000 0.1584699 0.171123
## 440649 0.1584699 1.0000000 0.646789
## 439161 0.1711230 0.6467890 1.000000
```

#### Convert adjacency matrix to an edge list

```r
#extracted connections for an undirected an graph
adjacency_edgeList(res)
```

```
##   source target             value
## 4     51 440649 0.158469945355191
## 7     51 439161 0.171122994652406
## 8 440649 439161 0.646788990825688
```

```r
#full matrix
adjacency_edgeList(res,symmetric=FALSE,diagonal=TRUE)
```

```
##   source target             value
## 1     51     51                 1
## 2 440649     51 0.158469945355191
## 3 439161     51 0.171122994652406
## 4     51 440649 0.158469945355191
## 5 440649 440649                 1
## 6 439161 440649 0.646788990825688
## 7     51 439161 0.171122994652406
## 8 440649 439161 0.646788990825688
## 9 439161 439161                 1
```
