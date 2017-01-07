#' @title CID_SDF
#' @param cid list of PubChem identifiers
#' @import dplyr ChemmineR 
CID_SDF<-function(cids,query.limit=25,...){
  
  #retrieve metabolite SDF from DB
  #DB should be a list with cids as names
  #for all missing in DB, look up using PubChem PUG
  #if update then update DB with cid entries
  #return list of SDF files for each cid
  
  # make sure all are numeric
  obj<-as.numeric(as.character(unlist(cids))) %>%
    .[!duplicated(.)] %>% na.omit()

   if(length(obj) == 0)  stop("Please supply valid input as numeric PubChem CIDs")
  
  #due to url string size limit query query.limit sdf obj at a time
  blocks<-as.list(split(obj, ceiling(seq_along(1:length(obj))/query.limit)))
  SDF<-lapply(1:length(blocks),function(i){
    url<-paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste(blocks[[i]] %>% unlist(),collapse=","),"/SDF")
    #format for output
    tmp<-read_sdf(url) %>% unclass(.)
    names(tmp)<-blocks[[i]]
    return(tmp)
  }) 
  
  #need to flatten list and name as cid from: http://stackoverflow.com/questions/19734412/flatten-nested-list-into-1-deep-list
  flatlist <- function(mylist){
    lapply(rapply(mylist, enquote, how="unlist"), eval)
  }
  
  #
 flatlist(SDF)
   
}  

#' @title read_sdf
read_sdf<-function (sdfstr) {
  
  #number of queries controlled in url
  if (length(sdfstr) > 1) {
    mysdf <- sdfstr
  } else {
    mysdf <- readLines(sdfstr)
  }
  
  y <- regexpr("^\\${4,4}", mysdf, perl = TRUE)
  index <- which(y != -1)
  indexDF <- data.frame(start = c(1, index[-length(index)] + 
                                    1), end = index)
  mysdf_list <- lapply(seq(along = indexDF[, 1]), function(x) mysdf[seq(indexDF[x, 
                                                                                1], indexDF[x, 2])])
  if (class(mysdf_list) != "list") {
    mysdf_list <- list(as.vector(mysdf_list))
  }
  names(mysdf_list) <- 1:length(mysdf_list)
  #mysdf_list <- new("SDFstr", a = mysdf_list)
  return(mysdf_list)
}

#' @import ChemmineR
SDF_tanimoto<-function(cmpd.DB){  
  #convert to SDFstr
  #depends on ChemmineR 
  require(ChemmineR)
  cmpd.sdf.list<-new("SDFstr", a = cmpd.DB)
  sd.list<-as(cmpd.sdf.list, "SDFset")
  cid(sd.list) <- sdfid(sd.list)
  
  # Convert base 64 encoded fingerprints to character vector, matrix or FPset object
  fpset <- fp2bit(sd.list, type=2)

  out<-sapply(rownames(fpset), function(x) ChemmineR::fpSim(x=fpset[x,], fpset,sorted=FALSE)) 
  obj<-as.matrix(out)
  
  return(obj)
}

#' @title CID_tanimoto
#' @param cid PubChem CID
#' @import ChemmineR dplyr
#' @details ... passed to CID_SDF
#' @export
CID_tanimoto<-function(cids,as=c("adjacency","edge.list"),...){
  cmpd.DB<-CID_SDF(cids,...)
  res<-SDF_tanimoto(cmpd.DB)
  if(as == 'edge.list') adjacency_edgeList(res) else res
} 

#' @title CID_tanimoto
#' @param mat adjacency matrix
#' @param symmetric TRUE if undirected keep lower.tri
#' @param diagonal TRUE to keep self-connections
#' @param mat adjacency matrix
#' @import reshape2
#' @details convert adjacency matrix into an edge list
#' @export
adjacency_edgeList<-function(mat,symmetric=TRUE,diagonal=FALSE){
  # bit of a hack around handling NA
  mat<-as.matrix(mat)
  id<-is.na(mat) # used to allow missing
  mat[id]<-"nna"
  if(symmetric){mat[lower.tri(mat)]<-"na"} # use to allow missing values
  if(!diagonal){diag(mat)<-"na"}
  obj<-melt(mat)
  colnames(obj)<-c("source","target","value")
  obj<-obj[!obj$value=="na",]
  obj$value[obj$value=="nna"]<-NA
  return(obj)
}	

#' @title get_KEGG_pairs
#' @param type of return 'main' type of reactions or 'full' everything
#' @param file RPAIRS list defaults to package data
#' @export
get_KEGG_pairs<-function(type="main",file=system.file("data/KEGG_RPAIRS", package = "metabomapr")){ 
  
  full<-read.table(file,sep="\t")
  
  if(type =="main"){
    out<-full[agrep("main",full[,3]),1:2] # now will get all main types
  } 

  if(type =="full"){
    out<-full
  }  
  return(as.matrix(out))
}

# #filter KEGG connections based on input
# input<-c("C15973","C00026","C05381","C15972","C00091","C00042")
# el<-get_KEGG_pairs()  
# match_edgeList(input,el)

#' @title get_KEGG_edgeList
#' @param input vector of kegg ids
#' @param el edge list of KEGG RPAIRs
#' @import dplyr
#' @export
get_KEGG_edgeList<-function(input,el=get_KEGG_pairs()){
  id<-as.character(el[,1]) %in% input  | as.character(el[,2]) %in% input
  tmp<-el[id,,drop=FALSE]
  res<-lapply(seq_along(input), function(i){
    id<-tmp[,1] %in% input[i] | tmp[,2] %in% input[i]
    tmp2<-tmp[id,,drop=FALSE]
    # lapply(tmp2, function(x) {
      tmp2[tmp2[,1] %in% input[-i] | tmp2[,2] %in% input[-i],]
    # }) %>% do.call("rbind",.)
  }) %>% do.call("rbind",.)
  rownames(res)<-NULL
  #TODO prevent duplicates above
  id<-apply(res,1,paste, collapse="_")
  res[!duplicated(id),] %>%
    data.frame(.,stringsAsFactors = FALSE) %>%
    setNames(.,c('source','target'))
  #doesn't take care of superimposed
}


#create shared edgelist
# db<-data.frame(ID=1:2,KEGG=c('C01416','C12448'),CID=c(4603,582838))
# end<-'ID'
# start<-'KEGG'
# edge_list<-data.frame(source=c('C12448','C01416'),target=c('C01416','C12448'),weight=1)
#convert_edgeIndex(edge_list,start,end,db)

#' @title convert_edgeIndex
#' @param edge_list matrix or data frame columns > 2 ignored
#' @param start column name of index matching edge list
#' @param end column name of index matching db
#' @param db data.frame containing mapping between start and end
#' @return edge list converted from start to end 
#' @import dplyr
#' @export 
convert_edgeIndex<-function(edge_list,start,end,db){
    
    #convert source/target to character for join  
    edge_list$source<-as.character(edge_list$source)
    edge_list$target<-as.character(edge_list$target)
    
    s<-left_join(edge_list[,'source',drop=FALSE],db,by=c('source' = start)) %>%
      select(one_of(end)) %>%
      setNames(.,'source')
    t<-left_join(edge_list[,'target',drop=FALSE],db,by=c('target' = start)) %>%
      select(one_of(end)) %>%
      setNames(.,'target')
  
    data.frame(s,t) 
} 
  
#' @title test_data
#' @return biochemical demo data
#' @export 
test_data<-function(){
  TCA.kegg <- c("C15973","C00026","C05381","C15972","C00091","C00042","C05379","C00311","C00036","C00024","C00149","C00417","C00158","C00022","C05125","C16254","C00122","C16255","C00074")
  TCA.names<-c('Enzyme N6-(dihydrolipoyl)lysine',	'2-Oxoglutarate',	'3-Carboxy-1-hydroxypropyl-ThPP',	'Enzyme N6-(lipoyl)lysine',	'Succinyl-CoA',	'Succinate',	'Oxalosuccinate',	'Isocitrate',	'Oxaloacetate',	'Acetyl-CoA',	'(S)-Malate',	'cis-Aconitate',	'Citrate',	'Pyruvate',	'2-(alpha-Hydroxyethyl)thiamine diphosphate',	'[Dihydrolipoyllysine-residue succinyltransferase] S-succinyldihydrolipoyllysine',	'Fumarate',	'[Dihydrolipoyllysine-residue acetyltransferase] S-acetyldihydrolipoyllysine',	'Phosphoenolpyruvate')
  TCA.CID <- c("[]","51", "440649","[]", "439161",   "1110",    "972",      "1198",     "970",      "6302",     "222656",   "643757",  "19782904", "1060",     "440568"  ,"[]", "21883788" ,"[]","1005" )
  data.frame(id=paste0('N',1:length(TCA.kegg)), KEGG=TCA.kegg,CID=TCA.CID,name=TCA.names,stringsAsFactors=FALSE)
}


#dev and test  
tests<-function(){
  library("metabomapr")
  
  #make sure main data inputs are character
  main<-test_data() 
  
  #### Convert PubChem CID to 
  #tanimoto similarity adjacency list
  type<-'CID'
  id<-main[,type]
  el<-CID_tanimoto(id,as='edge.list')
  cid_el<- el %>%
    convert_edgeIndex(.,start=type,end='id',db=main) %>%
    data.frame(.,el %>% select(-source,-target),type=type)
  
  ### get KEGG based biochemical connections
  type<-'KEGG'
  id<-main[,type]
  el<-get_KEGG_edgeList(id)
  kegg_el<- el %>%
    convert_edgeIndex(.,start=type,end='id',db=main) %>%
    data.frame(.,el %>% select(-source,-target),value=1,type=type)
  
  #combined edge list
  el<-rbind(kegg_el,cid_el)
  
}  
  