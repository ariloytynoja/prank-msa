library(XML)
par(mfrow=c(3,1))
alignment = xmlTreeParse("CAV2.xml")
root = xmlRoot(alignment)
for(i in 1:xmlSize(root[["nodes"]])) {
  if(xmlName(root[["nodes"]][[i]])=="node") {
   probs = strsplit(xmlValue(root[["nodes"]] [[i]][[2]]),",")[[1]]
   plot(as.numeric(probs)/100,type="h")
  }
}
