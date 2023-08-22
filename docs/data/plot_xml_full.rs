library(XML)

postscript(file="cav2.ps",horizontal=T,width=8,height=1.75)
par(mfrow=c(3,1),mar=c(1, 2.5, 0.5, 0.5),cex.axis=0.55,cex.lab=0.6,cex.main=1.5,tcl=0,mgp=c(0.8,0,0),lend=2,xaxs="i",yaxs="i")

data = xmlTreeParse("CAV2.xml")
root = xmlRoot(data)

nnodes = xmlSize(root[["nodes"]])
for(i in 1:nnodes) {
  if(xmlName(root[["nodes"]][[i]])=="node") {
    scores = as.numeric(strsplit(xmlValue(root[["nodes"]][[i]][[2]]),",")[[1]])
    plot(scores/100,xlim=range(c(0,length(scores))),ylim=range(c(0,1.1)),type="h",col="gray70",xlab="",ylab="P(slow)")
    text(300,0.85,xmlAttrs(root[["nodes"]][[i]])[["id"]],cex=0.7)
    lines(c(2284,2434),c(1.025,1.025))
    lines(c(2825,3012),c(1.025,1.025))
    lines(c(9018,9165),c(1.025,1.025))
  }
}

dev.off()
