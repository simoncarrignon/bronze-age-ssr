KWPV=isotopes[ isotopes$Site.abreviation %in% c("NKH","NPW"),]
IDVC=c( "SEALIP/LA/PBL/1", "SEALIP/LA/PBL/42", "SEALIP/LA/PBL/2", "SEALIP/LA/PBL/6", "SEALIP/LA/PBL/5", "SEALIP/LA/PBL/4", "SEALIP/LA/PBL/7", "SEALIP/LA/PBL/34", "SEALIP/LA/PBL/32", "SEALIP/LA/PBL/12", "SEALIP/LA/PBL/35", "SEALIP/LA/PBL/13", "SEALIP/LA/PBL/40", "SEALIP/LA/PBL/11", "SEALIP/LA/PBL/39", "SEALIP/LA/PBL/36", "SEALIP/LA/PBL/18", "SEALIP/LA/PBL/14", "SEALIP/LA/PBL/37", "SEALIP/LA/PBL/33", "SEALIP/LA/PBL/3", "SEALIP/LA/PBL/16", "SEALIP/LA/PBL/15", "SEALIP/LA/PBL/20", "SEALIP/LA/PBL/19", "SEALIP/LA/PBL/21", "SEALIP/LA/PBL/10", "SEALIP/LA/PBL/47")
VC=isotopes[ isotopes$LABEL.analytical %in% IDVC,]
sitecol=adjustcolor(RColorBrewer::brewer.pal(ncol,"Set1"),.8)


208Pb/204Pb:0.00365285642535325
207Pb/204Pb:0.00213129329946343
206Pb/204Pb:0.0160268856307887
cat(paste(ratioed(names(threshold)),paste0(round(as.numeric(threshold)*100,digit=3),"%"),sep=":"),sep="\n")



#FIGURE 6
pdf("Figure6.pdf",width=16,height=8)

par(mfrow=c(1,2),cex=1.4,mar=c(5,4,1,1))
for(rs in list(c(3,2),c(3,1))){
plot(rbind(VC[,ratioiso[rs]],KWPV[,ratioiso[rs]]),type="n",xlab=ratioed(ratioiso[rs[1]]),ylab=ratioed(ratioiso[rs[2]]))
points(KWPV[,ratioiso[rs]],bg=sitecol[1],pch=21)
points(VC[,ratioiso[rs]],bg=sitecol[2],pch=21)
points(KWPV[KWPV$LABEL.analytical %in%  c("SEALIP/TH/NPW/11","SEALIP/TH/NPW/1"),ratioiso[rs]],col="red",cex=2.5,pch=4,lwd=4)
}
val=legend("center",fill=sitecol,legend=c("KWPV","VC"),title="Complex",bty="n",ncol=1,cex=1)

dev.off()

KWPV=KWPV[!(KWPV$LABEL.analytical %in%  c("SEALIP/TH/NPW/11","SEALIP/TH/NPW/1")),]




threshold  <- sapply(ratioiso,function(r){
                     kpwd=getAllDist(KWPV[,r])
                     VC=getAllDist(VC[,r])
                     allconsistent=c(kpwd,VC)
                     quantile(allconsistent,probs=.95)[[1]]
})
allconsistents  <- sapply(ratioiso,function(r){
                     kpwd=getAllDist(KWPV[,r])
                     VC=getAllDist(VC[,r])
                     allconsistent=c(kpwd,VC)
})

par(mfrow=c(1,3),oma=c(0,0,2,0))
for(r in ratioiso){
    alld=allconsistents[[r]]
    plot(density(alld,from=0),main="",ylab="",xlab=paste0("ratio",gsub("\\.","/",gsub("X","",r))))
    #abline(v=threshold[r],lty=2,col="red")
    abline(v=quantile(alld,probs=c(.75,.85,.95)),lty=3,col="blue")
    mtext(paste0(c(75,85,95),"%"),3,0,at=quantile(alld,probs=c(.75,.85,.95)),col="blue",cex=.8)
    #mtext("oli's initial values",3,0,at=threshold[r],col="red")
    mtext("Distribution of distances for KWPV for all three ratios",3,0,outer=T)
}
mtext("Distribution of distances for KWPV for all three ratio",3,0,outer=T)

allVCs  <- lapply(ratioiso,function(r) getAllDist(VC[,r]))
allKWPVs  <- lapply(ratioiso,function(r) getAllDist(KWPV[,r]))
par(mfrow=c(3,3))
for(i in 1:3){ plot(density(allKWPVs[[i]])); plot(density(allVCs[[i]])); plot(density(c(allVCs[[i]],allKWPVs[[i]])));}

## FIGURE 15

#png("Figure15.png",width=1600,height=800,pointsize=18)
pdf("Figure15.pdf",width=16,height=8)
par(mfrow=c(1,2),mar=c(5,4,1,1),xpd=T,oma=rep(0,4))
gpe=paste0("gpe",tmpiso.clus$membership)
cnt=table(gpe)
ncol=length(cnt[cnt>1])
sitecol=adjustcolor(RColorBrewer::brewer.pal(ncol,"Set2"),.8)
names(sitecol)=names(cnt[cnt>1])
if(length(cnt)>sum(cnt>1)){sitecol=c(sitecol,"white");names(sitecol)[length(sitecol)]="No Consistency"}

plot(tmpiso[,ratioiso[3:2]],xlab=ratioed(ratioiso[3]),ylab=ratioed(ratioiso[2]),bg=sitecol[gpe],pch=21,cex=2,type="n")
text(tmpiso[,ratioiso[c(3,2)]],tmpiso$LABELS.shorts,cex=.4,pos=2)
points(tmpiso[,ratioiso[3:2]],bg=sitecol[gpe],pch=21,cex=2)
plot(tmpiso[,ratioiso[c(3,1)]],xlab=ratioed(ratioiso[3]),ylab=ratioed(ratioiso[1]),bg=sitecol[gpe],pch=21,cex=2,type="n")
text(tmpiso[,ratioiso[c(3,1)]],tmpiso$LABELS.shorts,,cex=.4,pos=2)
points(tmpiso[,ratioiso[c(3,1)]],bg=sitecol[gpe],pch=21,cex=2)
legend("topleft",legend=names(sitecol),pt.bg=sitecol,pch=21,pt.cex=2,ncol=2,cex=1)
dev.off()

#different order

par(mfrow=c(1,2),mar=c(0,0,0,1),oma=c(5,4,1,0),xpd=NA)
gpe=paste0("gpe",tmpiso.clus$membership)
cnt=table(gpe)
ncol=length(cnt[cnt>1])
sitecol=adjustcolor(RColorBrewer::brewer.pal(ncol,"Set2"),.8)
names(sitecol)=names(cnt[cnt>1])
if(length(cnt)>sum(cnt>1)){sitecol=c(sitecol,"white");names(sitecol)[length(sitecol)]="No Consistency"}

plot(tmpiso[,ratioiso[2:3]],xlab=ratioed(ratioiso[2]),ylab=ratioed(ratioiso[2]),bg=sitecol[gpe],pch=21,cex=2,type="n")
text(tmpiso[,ratioiso[c(2,3)]],tmpiso$LABELS.shorts,cex=.4,pos=4)
points(tmpiso[,ratioiso[2:3]],bg=sitecol[gpe],pch=21,cex=2)
plot(tmpiso[,ratioiso[c(1,3)]],xlab=ratioed(ratioiso[1]),ylab="",bg=sitecol[gpe],pch=21,cex=2,type="n",yaxt="n")
text(tmpiso[,ratioiso[c(1,3)]],tmpiso$LABELS.shorts,,cex=.4,pos=4)
points(tmpiso[,ratioiso[c(1,3)]],bg=sitecol[gpe],pch=21,cex=2)
legend("bottomright",legend=names(sitecol),pt.bg=sitecol,pch=21,pt.cex=2,ncol=2,cex=1)
###############


#FIGURE 16
tmpiso=isotopes[grep("Bronze",isotopes$Large.period),]

tmpiso.consmat=createRatioMatrix(tmpiso,ratioiso,threshold) 
dimnames(tmpiso.consmat)=list(tmpiso$LABELS.shorts,tmpiso$LABELS.shorts)
tmpiso.consmat[tmpiso.consmat<3]=0
tmpiso.consmat[tmpiso.consmat==3]=1
tmpiso.graph = graphFromAdj(tmpiso.consmat,label=tmpiso$LABELS.shorts)
tmpiso.clus=cluster_leiden(tmpiso.graph,objective_function="modularity")
tmpiso.sitemat=createAreaAdjMat(tmpiso,tmpiso$Site.abreviation,tmpiso.clus$membership,match=RGmatchNormalised)
tmpiso.sitegraph = graphFromAdj(tmpiso.sitemat,label=rownames(tmpiso.sitemat))

tmpiso.sitecluster=cluster_leiden(tmpiso.sitegraph,objective_function = "modularity")
V(tmpiso.sitegraph)$color=tmpiso.sitecluster$membership
V(tmpiso.sitegraph)$size=20
#V(tmpiso.sitegraph)$size=betweenness(tmpiso.sitegraph)
E(tmpiso.sitegraph)$width=2*(exp(E(tmpiso.sitegraph)$weight)-1)
tmpiso.sitegraph=colorEdges(tmpiso.sitegraph)
pdf("Figure16.pdf",width=8,height=8)
par(mar=rep(0,4))
plot(tmpiso.sitegraph,main="")
dev.off()
###############



#Figure 17
V(tmpiso.sitegraph)$label.dist=5
tmpiso.sitegraph.coords=coordForGraph(tmpiso.sitegraph,tmpiso,"Site.abreviation")
invalidcoords=apply(apply(tmpiso.sitegraph.coords,2,is.na),1,any)
tmpiso.sf=st_as_sf(tmpiso[!is.na(tmpiso[,7]),],coords=6:7,crs=st_crs(sealip))
fullgraph.sp=st_as_sf(as.data.frame(tmpiso.sitegraph.coords[!invalidcoords,]),coords=1:2,crs=st_crs(sealip))
netssr.sf=st_as_sf(as.data.frame(coords[!is.na(coords[,1]),]),coords=c(1,2))
lims=st_bbox(st_buffer(st_geometry(tmpiso.sf),100000))

pdf("Figure17.pdf",width=9,height=8)
par(mar=rep(0,4))
plot(st_geometry(st_crop(sealip,st_buffer(netssr.sf,50000))),col=adjustcolor("pink",.2),extent=lims)
E(tmpiso.sitegraph)$color[E(tmpiso.sitegraph)$color=="black"]=adjustcolor("black",.2)
plot(tmpiso.sitegraph,layout=tmpiso.sitegraph.coords,add=T,rescale=F)
plot(netssr,layout=coords[V(netssr)$name,],add=T,rescale=F)
box(which = "outer")
dev.off()
        invalidcoords=apply(apply(tmpiso.sitegraph.coords,2,is.na),1,any)

pdf(paste0("Figure18.pdf"),width=21,height=8)
par(mfrow=c(1,3),mar=rep(0,4),oma=rep(0,4))
for(gp in which(lengths(tmpiso.sitecluster,use.names=F)>1)){
    sub=subgraph(tmpiso.sitegraph,tmpiso.sitecluster[[gp]])
    sub=set_graph_attr(sub,"layout",tmpiso.sitegraph.coords[V(sub)$name,])
    #pdf(paste0("Figure17_",gp,".pdf"),width=9,height=8)
    #try({
    #E(sub)$width=3
    graph.sp=st_as_sf(as.data.frame(tmpiso.sitegraph.coords[!invalidcoords,]),coords=1:2,crs=st_crs(sealip))
    plot(st_geometry(st_crop(sealip,st_buffer(netssr.sf,50000))),col=adjustcolor("pink",.2),extent=lims)
    plot(sub,add=T,rescale=F)
    plot(netssr,layout=coords[V(netssr)$name,],add=T,rescale=F)
    box(which="figure" ) 
    #})
    #dev.off()
}
dev.off()


