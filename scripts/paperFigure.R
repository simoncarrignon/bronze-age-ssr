source("R/tools.R")
isotopes=read.csv("data/DATA isotopie_Pryce_V4.2.1.csv")
KWPV=isotopes[ isotopes$Site.abreviation %in% c("NKH","NPW"),]
IDVC=c( "SEALIP/LA/PBL/1", "SEALIP/LA/PBL/42", "SEALIP/LA/PBL/2", "SEALIP/LA/PBL/6", "SEALIP/LA/PBL/5", "SEALIP/LA/PBL/4", "SEALIP/LA/PBL/7", "SEALIP/LA/PBL/34", "SEALIP/LA/PBL/32", "SEALIP/LA/PBL/12", "SEALIP/LA/PBL/35", "SEALIP/LA/PBL/13", "SEALIP/LA/PBL/40", "SEALIP/LA/PBL/11", "SEALIP/LA/PBL/39", "SEALIP/LA/PBL/36", "SEALIP/LA/PBL/18", "SEALIP/LA/PBL/14", "SEALIP/LA/PBL/37", "SEALIP/LA/PBL/33", "SEALIP/LA/PBL/3", "SEALIP/LA/PBL/16", "SEALIP/LA/PBL/15", "SEALIP/LA/PBL/20", "SEALIP/LA/PBL/19", "SEALIP/LA/PBL/21", "SEALIP/LA/PBL/10", "SEALIP/LA/PBL/47")
VC=isotopes[ isotopes$LABEL.analytical %in% IDVC,]


#208Pb/204Pb:0.00365285642535325
#207Pb/204Pb:0.00213129329946343
#206Pb/204Pb:0.0160268856307887

ratioiso=c("X208Pb.204Pb","X207Pb.204Pb","X206Pb.204Pb")
threshold=readRDS("data/thresholds.RDS")
cat(paste(ratioed(names(threshold)),paste0(round(as.numeric(threshold)*100,digit=3),"%"),sep=":"),sep="\n")



#FIGURE 6
pdf("Figure6.pdf",width=16,height=8)
sitecol=1:2
par(mfrow=c(1,2),cex=1.4,mar=c(5,4,1,1))
for(rs in list(c(3,2),c(3,1))){
plot(rbind(VC[,ratioiso[rs]],KWPV[,ratioiso[rs]]),type="n",xlab=ratioed(ratioiso[rs[1]]),ylab=ratioed(ratioiso[rs[2]]))
points(KWPV[,ratioiso[rs]],bg=sitecol[1],pch=21)
points(VC[,ratioiso[rs]],bg=sitecol[2],pch=21)
points(KWPV[KWPV$LABEL.analytical %in%  c("SEALIP/TH/NPW/11","SEALIP/TH/NPW/1"),ratioiso[rs]],col="red",cex=2.5,pch=4,lwd=4)
}
val=legend("center",fill=sitecol,legend=c("KWPV","VC"),title="Complex",bty="n",ncol=1,cex=1)

dev.off()

## FIGURE 15
tmpiso=isotopes[grep("Bronze",isotopes$Large.period..SEA.equivalent.),]
library(igraph)
tmpiso.sitegraph=readRDS(file="data/tmpiso.sitegraph.RDS")
tmpiso.graph=readRDS(file="data/tmpiso.graph.RDS")
tmpiso.clus=readRDS(file="data/tmpiso.clus.RDS")

#png("Figure15.png",width=1600,height=800,pointsize=18)
pdf("Figure15.pdf",width=16,height=8)
par(mfrow=c(1,2),mar=c(5,4,1,1),xpd=T,oma=rep(0,4))
gpe=paste0("gpe",tmpiso.clus$membership)
cnt=table(gpe)
ncol=length(cnt[cnt>1])
sitecol=adjustcolor(RColorBrewer::brewer.pal(ncol,"Set2"),.8)
pchs=rep(1,length(sitecol))
if(length(sitecol)<length(cnt[cnt>1])){
          pchs[(length(sitecol)+1):length(cnt[cnt>1])]=2
          sitecol=rep(sitecol,2)[1:length(cnt[cnt>1])]
}
names(sitecol)=names(cnt[cnt>1])
names(pchs)=names(sitecol)
pchs=c(pchs,rep(1,length(cnt[cnt<=1])))
names(pchs)[names(pchs)==""]=names(cnt[cnt<=1])
#pchs["gpe1"]=2
#pchs["gpe6"]=1

if(length(cnt)>sum(cnt>1)){sitecol=c(sitecol,"white");names(sitecol)[length(sitecol)]="No Consistency";}

plot(tmpiso[,ratioiso[3:2]],xlab=ratioed(ratioiso[3]),ylab=ratioed(ratioiso[2]),type="n")
text(tmpiso[,ratioiso[c(3,2)]],tmpiso$LABELS.shorts,cex=.4,pos=2)
points(tmpiso[,ratioiso[3:2]],bg=sitecol[gpe],pch=20+pchs[gpe],cex=2)
plot(tmpiso[,ratioiso[c(3,1)]],xlab=ratioed(ratioiso[3]),ylab=ratioed(ratioiso[1]),type="n")
text(tmpiso[,ratioiso[c(3,1)]],tmpiso$LABELS.shorts,,cex=.4,pos=2)
points(tmpiso[,ratioiso[c(3,1)]],bg=sitecol[gpe],pch=20+pchs[gpe],cex=2)
legend("topleft",legend=names(sitecol),pt.bg=sitecol,pt.cex=2,ncol=2,cex=1,pch=pchs+20)
dev.off()



#FIGURE 16
pdf("Figure16.pdf",width=8,height=8)
par(mar=rep(0,4))
plot(tmpiso.sitegraph,main="")
dev.off()




#Figure 17
V(tmpiso.sitegraph)$label.dist=5
tmpiso.sitegraph.coords=coordForGraph(tmpiso.sitegraph,tmpiso,"Site.abreviation")
invalidcoords=apply(apply(tmpiso.sitegraph.coords,2,is.na),1,any)
library(sf)
sealip=st_read("maps/allseaborder_lowq.shp")
tmpiso.sf=st_as_sf(tmpiso[!is.na(tmpiso[,7]) | !is.na(tmpiso[,8]),],coords=7:8,crs=st_crs(sealip))
fullgraph.sp=st_as_sf(as.data.frame(tmpiso.sitegraph.coords[!invalidcoords,]),coords=1:2,crs=st_crs(sealip))
coords=readRDS(file="data/coordnetssr.RDS")
netssr=readRDS(file="data/netssr.RD")
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

tmpiso.sitecluster=readRDS("data/tmpiso.sitecluster.RDS")

pdf(paste0("Figure18.pdf"),width=21,height=8)
par(mfrow=c(1,4),mar=rep(0,4),oma=rep(0,4))
for(gp in which(lengths(tmpiso.sitecluster,use.names=F)>1)){
    sub=subgraph(tmpiso.sitegraph,tmpiso.sitecluster[[gp]])
    sub=set_graph_attr(sub,"layout",tmpiso.sitegraph.coords[V(sub)$name,])
    try({
    E(sub)$width=3
    graph.sp=st_as_sf(as.data.frame(tmpiso.sitegraph.coords[!invalidcoords,]),coords=1:2,crs=st_crs(sealip))
    plot(st_geometry(st_crop(sealip,st_buffer(netssr.sf,50000))),col=adjustcolor("pink",.2),extent=lims)
    plot(sub,add=T,rescale=F)
    plot(netssr,layout=coords[V(netssr)$name,],add=T,rescale=F)
    box(which="figure" ) 
    })
}
dev.off()

for(gp in which(lengths(tmpiso.sitecluster,use.names=F)>1)){
    sub=subgraph(tmpiso.sitegraph,tmpiso.sitecluster[[gp]])
    sub=set_graph_attr(sub,"layout",tmpiso.sitegraph.coords[V(sub)$name,])
    pdf(paste0("Figure17_",gp,".pdf"),width=9,height=8)
    par(mar=rep(0,4),oma=rep(0,4))
    try({
    E(sub)$width=3
    graph.sp=st_as_sf(as.data.frame(tmpiso.sitegraph.coords[!invalidcoords,]),coords=1:2,crs=st_crs(sealip))
    plot(st_geometry(st_crop(sealip,st_buffer(netssr.sf,50000))),col=adjustcolor("pink",.2),extent=lims)
    plot(sub,add=T,rescale=F)
    plot(netssr,layout=coords[V(netssr)$name,],add=T,rescale=F)
    box(which="figure" ) 
    })
    dev.off()
}

## Script to convert in png high quality:
#for i in {1..4} ; do convert -density 600 Figure17_$i.pdf Figure17_$i.png ; done
#for i in {15..18} ; do convert -density 600 Figure$i.pdf Figure$i.png ; done
