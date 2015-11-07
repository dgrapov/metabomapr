
### Install

```r
#devtools::install_github("dgrapov/metabomapr")
library("metabomapr")
```

### Examples
#### Convert PubChem CID to tanimoto similarity adjacency matrix

```r
cids<-c("[]","51", "440649","[]", "439161",  NA )
CID_tanimoto(cids)
```

```
##               51    440649   439161
## 51     1.0000000 0.1584699 0.171123
## 440649 0.1584699 1.0000000 0.646789
## 439161 0.1711230 0.6467890 1.000000
```
