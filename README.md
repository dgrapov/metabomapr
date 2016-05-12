
### Metabolomic network edge list generation. Currently includes:
* structural similarity
* KEGG biochemical connections

#### See [MetaMapR](http://dgrapov.github.io/MetaMapR/) for more features.

### Install

```r
#devtools::install_github("dgrapov/metabomapr")
library("metabomapr")

TCA.kegg <- c("C15973","C00026","C05381","C15972","C00091","C00042","C05379","C00311","C00036","C00024","C00149","C00417","C00158","C00022","C05125","C16254","C00122","C16255","C00074")
TCA.names<-c('Enzyme N6-(dihydrolipoyl)lysine',	'2-Oxoglutarate',	'3-Carboxy-1-hydroxypropyl-ThPP',	'Enzyme N6-(lipoyl)lysine',	'Succinyl-CoA',	'Succinate',	'Oxalosuccinate',	'Isocitrate',	'Oxaloacetate',	'Acetyl-CoA',	'(S)-Malate',	'cis-Aconitate',	'Citrate',	'Pyruvate',	'2-(alpha-Hydroxyethyl)thiamine diphosphate',	'[Dihydrolipoyllysine-residue succinyltransferase] S-succinyldihydrolipoyllysine',	'Fumarate',	'[Dihydrolipoyllysine-residue acetyltransferase] S-acetyldihydrolipoyllysine',	'Phosphoenolpyruvate')
TCA.CID <- c("[]","51", "440649","[]", "439161",   "1110",    "972",      "1198",     "970",      "6302",     "222656",   "643757",  "19782904", "1060",     "440568"  ,"[]", "21883788" ,"[]","1005" )
main<-data.frame(id=1:length(TCA.kegg), KEGG=TCA.kegg,CID=TCA.CID,name=TCA.names)
```

### Examples
#### Convert PubChem CID to tanimoto similarity adjacency list

```r
cids<-main$CID
res1<-CID_tanimoto(cids,as='edge.list')
head(res1)
```

```
##    source target             value
## 16     51 440649 0.158469945355191
## 31     51 439161 0.171122994652406
## 32 440649 439161 0.646788990825688
## 46     51   1110 0.702702702702703
## 47 440649   1110 0.146892655367232
## 48 439161   1110 0.147540983606557
```

### get KEGG based biochemical connections

```r
kegg<-main$KEGG
res2<-get_KEGG_edgeList(kegg)
head(res2)
```

```
##      V1       V2      
## [1,] "C15972" "C15973"
## [2,] "C15973" "C16254"
## [3,] "C15973" "C16255"
## [4,] "C00026" "C00091"
## [5,] "C00026" "C00311"
## [6,] "C00026" "C05381"
```

### combine KEGG and CID edge lists using a common index

```r
library(dplyr)
#main db
head(main)
```

```
##   id   KEGG    CID                            name
## 1  1 C15973     [] Enzyme N6-(dihydrolipoyl)lysine
## 2  2 C00026     51                  2-Oxoglutarate
## 3  3 C05381 440649  3-Carboxy-1-hydroxypropyl-ThPP
## 4  4 C15972     []        Enzyme N6-(lipoyl)lysine
## 5  5 C00091 439161                    Succinyl-CoA
## 6  6 C00042   1110                       Succinate
```

```r
main$CID<-main$CID %>% as.character() %>% as.numeric()
cid_el<-res1 %>% 
  convert_edgeIndex(.,start='CID',end='id',db=main) %>%
  mutate(type="CID")

head(cid_el)
```

```
##   source target             value type
## 1      2      3 0.158469945355191  CID
## 2      2      5 0.171122994652406  CID
## 3      3      5 0.646788990825688  CID
## 4      2      6 0.702702702702703  CID
## 5      3      6 0.146892655367232  CID
## 6      5      6 0.147540983606557  CID
```

```r
kegg_el<-res2 %>% 
  convert_edgeIndex(.,start='KEGG',end='id',db=main) %>%
  mutate(value=1,type="KEGG")

head(kegg_el)
```

```
##   source target value type
## 1      4      1     1 KEGG
## 2      1     16     1 KEGG
## 3      1     18     1 KEGG
## 4      2      5     1 KEGG
## 5      2      8     1 KEGG
## 6      2      3     1 KEGG
```

```r
final_edgeList<-rbind(kegg_el,cid_el)
```

