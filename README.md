---
output:
  html_document:
    keep_md: yes
---

> Metabolomic network edge list generation

### Supported connections:
* structural similarity based on [pubChem CID](https://pubchem.ncbi.nlm.nih.gov/)
* biochemical connections based on [KEGG](http://www.genome.jp/kegg/)

#### See [MetaMapR](http://dgrapov.github.io/MetaMapR/) for more features.

### Examples

```r
#devtools::install_github("dgrapov/metabomapr")
library(metabomapr)
library(dplyr)
```

### Demo data

```r
main<-test_data() 
head(main)
```

```
##   id   KEGG    CID                            name
## 1 N1 C15973     [] Enzyme N6-(dihydrolipoyl)lysine
## 2 N2 C00026     51                  2-Oxoglutarate
## 3 N3 C05381 440649  3-Carboxy-1-hydroxypropyl-ThPP
## 4 N4 C15972     []        Enzyme N6-(lipoyl)lysine
## 5 N5 C00091 439161                    Succinyl-CoA
## 6 N6 C00042   1110                       Succinate
```

### Convert PubChem CID to tanimoto based chemical similarity

```r
type<-'CID'
id<-main[,type]
cid_el<-CID_tanimoto(id,as='edge.list')
head(cid_el)
```

```
##    source target     value
## 16     51 440649 0.1584699
## 31     51 439161 0.1711230
## 32 440649 439161 0.6467890
## 46     51   1110 0.7027027
## 47 440649   1110 0.1468927
## 48 439161   1110 0.1475410
```

### get KEGG based biochemical connections

```r
type<-'KEGG'
id<-main[,type]
kegg_el<-get_KEGG_edgeList(id)
head(kegg_el)
```

```
##   source target
## 1 C15972 C15973
## 2 C15973 C16254
## 3 C15973 C16255
## 4 C00026 C00091
## 5 C00026 C00311
## 6 C00026 C05381
```

### Combine edge lists using a common index

```r
type<-'CID'
cid_el<- cid_el %>%
  convert_edgeIndex(.,start=type,end='id',db=main) %>%
  data.frame(.,cid_el %>% select(-source,-target),type=type)

type<-'KEGG'
kegg_el<- kegg_el %>%
  convert_edgeIndex(.,start=type,end='id',db=main) %>%
  data.frame(.,kegg_el %>% select(-source,-target),value=1,type=type)

#combined edge list
el<-rbind(kegg_el,cid_el)
head(el)
```

```
##   source target value type
## 1     N4     N1     1 KEGG
## 2     N1    N16     1 KEGG
## 3     N1    N18     1 KEGG
## 4     N2     N5     1 KEGG
## 5     N2     N8     1 KEGG
## 6     N2     N3     1 KEGG
```

## Using opencpu API

Define server
```
options('open_cpu_url' = 'http://localhost/ocpu/')
options('metabomapr_url' = 'library/metabomapr/R/' )
```

Make call using simple [opencpuclient](https://github.com/dgrapov/ocpuclient)
```
library(ocpuclient)
type<-'CID'
id<-main[,type]

body<-list(type=type,cids=id,as='edge.list')

fun<-'CID_tanimoto'

x<-ocpuclient:::ocpu_post(fun,body=body,pkg_url = getOption('metabomapr_url'), base_url=getOption('open_cpu_url'))

x$results 
```


### About
* contact: createdatasol@gmail.com
* updated: 2021-09-21
