fun_range <- function(x,n=100) {                              # Create user-defined function
  (x - min(x)) / (max(x) - min(x)) * 100
}
simpmatch  <-  function(i,j){
    if(length(i)>length(j)){
        tmp=i
        i=j
        j=tmp
    }
    length(na.omit(match(i,j)))/(length(c(i,j)))
}

ratioSpan  <-  function(a,b=NULL){
    maxr=NULL
    minr=NULL
    if(length(a)>=2){
        maxr=max(a)
        minr=min(a)
    }
    else{
        maxr=max(a,b)
        minr=min(a,b)
    }
    minr/maxr
}

same <- function(a,b=NULL,threshold=0.01){
    maxr=minr=NULL
    if(length(a)==2){
        maxr=max(a)
        minr=min(a)
    }
    else{
        maxr=max(a,b)
        minr=min(a,b)
    }
    if(1-(minr/maxr) <= threshold) return(1)
    else return(0)
}

#cleanlist: dataframe
# ratioiso the column used to store isotope ratio
#threshold the threshold for each ratio
createRatioMatrix <- function(dataset,ratioiso,thresholds){
    adjmat=matrix(0,nrow=nrow(dataset),ncol=nrow(dataset))
    for(i in 1:(nrow(dataset)-1)){
        for(j in (i+1):nrow(dataset)){
            adjmat[j,i]=adjmat[i,j]=sum(sapply(ratioiso,function(iso)same(dataset[c(i,j),iso],threshold=thresholds[iso])))
        }
    }
    return(adjmat)
}
       
graphFromAdj <- function(adjmat,mask=NULL,label=NULL){
    require(igraph)
    if(!is.null(mask) &&  mask>0){
        adjmat[adjmat<mask]=0
        adjmat[adjmat>=mask]=1
    }
    graph=graph_from_adjacency_matrix(adjmat,diag=F,weighted=T,mode="undirected")
    V(graph)$label=label
    V(graph)$size=8
    V(graph)$frame.color="black"
    V(graph)$frame.width=.5
    V(graph)$label.color=adjustcolor("black",.6)
    E(graph)$color="black"
    E(graph)$width=.4
    #E(graph)$arrow.mode="-"
    graph=set_graph_attr(graph,"layout",layout_with_fr(graph))
    graph
}

#' @param dataset all data
#' @param sites how dataset should be grouped by site
#' @param match a function that coompute how sites match
creatSiteAdjMat <- function(dataset,sitesgroup,artgroup,match=NULL){
    sites=unique(sitesgroup)
    adjmat=matrix(NA,nrow=length(sites),ncol=length(sites))
    for(s1 in 1:length(sites)){
        for(s2 in 1:length(sites)){
            c1=artgroup[ sitesgroup == sites[s1]]
            c2=artgroup[ sitesgroup == sites[s2]]
            adjmat[s1,s2]=match(c1,c2)
        }
    }
    rownames(adjmat)=sites
    colnames(adjmat)=sites
    adjmat
}

groupGraphIso <- function(subset,coords,grouping,ratioiso,recomp){
    subset=subset[!is.na(subset[[coords[1]]]),]
    subset=subset[subset[[grouping]]!="",]
    group=unique(subset[,grouping])
    getSitesCoordiates=unique(subset[,c(grouping,coords)])
    coordsites=t(sapply(group,function(lbl)apply(unique(getSitesCoordiates[getSitesCoordiates[[grouping]] == lbl,coords]),2,mean)))
    allsites=createRatioMatrix(subset,ratioiso,recomp)
    graph=graphFromAdj(allsites,mask=3,label=subset[[grouping]])
    comus=cluster_louvain(graph)$memberships[1,]
    adjmat=creatSiteAdjMat(subset,subset[[grouping]],comus,match=simpmatch)
    groupgraph=graph.adjacency(adjmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    groupgraph=set_graph_attr(groupgraph,"layout",coordsites)
    groupcomu=cluster_louvain(groupgraph)$memberships[1,]
    print(length(unique(groupcomu)))
    E(groupgraph)$color=groupcomu[head_of(groupgraph,E(groupgraph))]
    V(groupgraph)$color=groupcomu
    V(groupgraph)$size=100
    E(groupgraph)$width=E(groupgraph)$weight*2
    V(groupgraph)$label.color="black"
    groupgraph
}

groupGraphElmt <- function(subset,coords,grouping,listelements){
    subset=subset[!is.na(subset[[coords[1]]]),]
    subset=subset[subset[[grouping]]!="",]
    group=unique(subset[,grouping])
    getSitesCoordiates=unique(subset[,c(grouping,coords)])
    coordsites=t(sapply(group,function(lbl)apply(unique(getSitesCoordiates[getSitesCoordiates[[grouping]] == lbl,coords]),2,mean)))
    elts=subset[,listelements]
    elts.pca=prcomp(elts)
    elts.dist=dist(elts.pca$x) #compute distance between PCA scores
    elts.ew=1/exp(as.matrix(elts.dist))
    elts.graph=graph.adjacency(elts.ew, mode = "undirected", weighted = TRUE, diag = FALSE)
    elts.commu=membership(cluster_louvain(elts.graph))
    comus=cluster_louvain(graph)$memberships[1,]
    adjmat=creatSiteAdjMat(subset,subset[[grouping]],elts.commu,match=simpmatch)
    groupgraph=graph.adjacency(adjmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    groupgraph=set_graph_attr(groupgraph,"layout",coordsites)
    groupcomu=cluster_louvain(groupgraph)$memberships[1,]
    print(length(unique(groupcomu)))
    E(groupgraph)$color=groupcomu[head_of(groupgraph,E(groupgraph))]
    V(groupgraph)$color=groupcomu
    V(groupgraph)$size=100
    E(groupgraph)$width=E(groupgraph)$weight*2
    V(groupgraph)$label.color="black"
    groupgraph
}
