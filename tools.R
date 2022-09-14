fun_range <- function(x,n=100) {                              # Create user-defined function
  (x - min(x)) / (max(x) - min(x)) * 100
}

# A normalised verison of the match used by Jelena and Miljana that cont how many similar elements are shared between two vectors
simpmatch  <-  function(i,j){
    if(length(i)>length(j)){
        tmp=i
        i=j
        j=tmp
    }
    length(na.omit(match(i,j)))/(length(c(i,j)))
}

divmatch  <-  function(i,j,normalize=18){
    if(length(i)>length(j)){
        tmp=i
        i=j
        j=tmp
    }
    length(na.omit(match(i,j)))
}

# return the ratio betwen the smallest and largest element of a and b 
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

# return tru if the ratio between a an b is beloew a threshold (overlap with ratioSpan)
same <- function(a,b=NULL,threshold=0.01,log=F){
    maxr=minr=NULL
    if(length(a)==2){
        maxr=max(a)
        minr=min(a)
    }
    else{
        maxr=max(a,b)
        minr=min(a,b)
    }
    if(log){
        minr=log(minr)
        maxr=log(maxr)
    }
    if(1-(minr/maxr) <= threshold) return(1)
    else return(0)
}

#cleanlist: dataframe
# ratioiso the column used to store isotope ratio
#threshold the threshold for each ratio
createRatioMatrix <- function(dataset,ratioiso,thresholds,log=F){
    adjmat=matrix(0,nrow=nrow(dataset),ncol=nrow(dataset))
    for(i in 1:(nrow(dataset)-1)){
        for(j in (i+1):nrow(dataset)){
            adjmat[j,i]=adjmat[i,j]=sum(sapply(ratioiso,function(iso)same(dataset[c(i,j),iso],threshold=thresholds[iso],log=log)))
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
#' @param area how dataset should be grouped by site
#' @param match a function that coompute how sites match
#' @param compoistion the value used to match the area
createAreaAdjMat <- function(dataset,area,composition,match=NULL){
    topLevelArea=unique(area)
    adjmat=matrix(NA,nrow=length(topLevelArea),ncol=length(topLevelArea))
    for(s1 in 1:length(topLevelArea)){
        for(s2 in 1:length(topLevelArea)){
            c1=composition[ area == topLevelArea[s1]]
            c2=composition[ area == topLevelArea[s2]]
            adjmat[s1,s2]=match(c1,c2)
        }
    }
    rownames(adjmat)=topLevelArea
    colnames(adjmat)=topLevelArea
    adjmat
}

#' @param subset all data
#' @param coords the colonms with coordinates in subset
#' @param recomp a list of the same size than `ratioiso` used to define the threhold below which two sample are similiar
#' @param ratioiso the list of isotops ratio to use
groupGraphIso <- function(subset,coords,grouping,ratioiso,recomp){
    subset=subset[!is.na(subset[[coords[1]]]),]
    subset=subset[subset[[grouping]]!="",]
    allsites=createRatioMatrix(subset,ratioiso,recomp)
    graph=graphFromAdj(allsites,mask=3,label=subset[[grouping]])
    comus=cluster_louvain(graph)$memberships[1,]
    groupedGraphAndCommunity(subset,coords,grouping,comus)
}

#' @param subset all data
#' @param coords the colonms with coordinates in subset
#' @param grouping the column use to define geographical sites
#' @param listelements the list of elements used to compute the distance between two samples
#' @return a graph
groupGraphElmt <- function(subset,coords,grouping,listelements){
    subset=subset[!is.na(subset[[coords[1]]]),]
    subset=subset[subset[[grouping]]!="",]
    elts=subset[,listelements]
    elts.pca=prcomp(elts)
    elts.dist=dist(elts.pca$x) #compute distance between PCA scores
    elts.ew=1/exp(as.matrix(elts.dist))
    elts.graph=graph.adjacency(elts.ew, mode = "undirected", weighted = TRUE, diag = FALSE)
    elts.commu=membership(cluster_louvain(elts.graph))
    groupedGraphAndCommunity(subset,coords,grouping,elts.commu)
}

#' @param dataset a `data.frame` storing all data
#' @param grouping the column in dataset used to group the data together
#' @param categories this is a liste of the size of unique(dataset$grouping) that tell how each group of grouping are linked
#' @param coords coordinates of the echantillon from original data `dataset`
groupedGraphAndCommunity <- function(dataset,coords,grouping,categories){
    group=unique(dataset[,grouping])
    getSitesCoordiates=unique(dataset[,c(grouping,coords)])
    coordsites=t(sapply(group,function(lbl)apply(unique(getSitesCoordiates[getSitesCoordiates[[grouping]] == lbl,coords]),2,mean)))
    adjmat=creatSiteAdjMat(dataset,dataset[[grouping]],categories,match=simpmatch)
    groupgraph=graph.adjacency(adjmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    groupgraph=set_graph_attr(groupgraph,"layout",coordsites)
    #should use graphFromAdj instead of the two previous lines
    groupcomu=cluster_louvain(groupgraph)$memberships[1,]
    E(groupgraph)$color=groupcomu[head_of(groupgraph,E(groupgraph))]
    V(groupgraph)$color=groupcomu
    V(groupgraph)$size=100
    E(groupgraph)$width=E(groupgraph)$weight*2
    V(groupgraph)$label.color="black"
    groupgraph
}


#'ultimately should be a plot that, on the usual biplot of isotopos ratio, should plot links stored in a consistency matrix (ie network adjacency matrix)
#' @param dataset a `data.frame` storing all data
#' @param categories a vector of the size of nrow(datase) that tell  group of each element
#' @param consistencmatrix wich points need to be linked
plotConsistensy <- function(dataset,col2plot,consistencmat,ratiolim=NULL,colour){
    par(mfrow=c(2,1))
    plot(dataset[,ratioiso[3]],dataset[,ratioiso[2]],xlim=ratiolim[[ratioiso[3]]],ylim=ratiolim[[ratioiso[2]]],xlab=ratioiso[3],ylab=ratioiso[2],bg=colour,pch=21)
    for(a in 1:nrow(dataset)){
        coorda=dataset[a,ratioiso[c(3,2)]]
        coordbs=dataset[consistencymat[a,]==1,ratioiso[c(3,2)],drop=F]
        if(nrow(coordbs)>1) apply(coordbs,1,function(b)segments(x0=coorda[[1]],y0=coorda[[2]],x1=b[[1]],y1=b[[2]],lwd=.1))
    }
    points(dataset[,ratioiso[3]],dataset[,ratioiso[2]],bg=colour,pch=21)
    text(dataset[,ratioiso[3]],dataset[,ratioiso[2]],dataset$Site.abreviation,cex=.7,pos=1)

    plot(dataset[,ratioiso[3]],dataset[,ratioiso[1]],xlim=ratiolim[[ratioiso[3]]],ylim=ratiolim[[ratioiso[1]]],xlab=ratioiso[3],ylab=ratioiso[1],bg=colour,pch=21)
    for(a in 1:nrow(dataset)){
        coorda=dataset[a,ratioiso[c(3,1)]]
        coordbs=dataset[consistencymat[a,]==1,ratioiso[c(3,1)],drop=F]
        if(nrow(coordbs)>1) apply(coordbs,1,function(b)segments(x0=coorda[[1]],y0=coorda[[2]],x1=b[[1]],y1=b[[2]],lwd=.1))
    }
    text(dataset[,ratioiso[3]],dataset[,ratioiso[1]],dataset$Site.abreviation,cex=.7,pos=1)
    points(dataset[,ratioiso[3]],dataset[,ratioiso[1]],bg=colour,pch=21)
    legend("bottomright",legend=levels(as.factor(dataset$Country)),pch=1:11)
}


#' a function to retrun coordinates mathcing grap[h label. If graph labels implies multiples coordinate sa mean will be done
#' only cthe coordinates are plotted and note a "layouted" graph to allow to play with the igraph graph and sf objet, without haveing to handle the NA (who prevent the creation of sf objects)  
#' @param categories a vector of the size of nrow(datase) that tell  group of each element
#' @param dataset a `data.frame` storing all data: need to containt, the coordinates associate to labels  \emph{and} the labels  that should match the graph labels.
#' @param labelscol: the corresponding colonm in the dataset where label are stored
#' @param coords: the corresponding colonm in the dataset where coordinates are stored
coordForGraph <- function(graph,dataset,labelscol,coords=c("Easting...Longitude","Northing.Latitude")){
    t(sapply(V(graph)$name,function(lbl)apply(unique(dataset[ dataset[[labelcol]] == lbl,coords]),2,mean,na.rm=T)))
    }
