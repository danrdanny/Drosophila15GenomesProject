#
library(ggplot2)

data <- read.table(sep='\t', header=TRUE, file="~/projects/GitHub/Drosophila15GenomesProject/dere.scaffold485.forR.tsv")

newRow1 <- data.frame(Scaffold="1",Start=0,End=1,Size=1,Asm="nanopore")
newRow2 <- data.frame(Scaffold="2",Start=0,End=1,Size=1,Asm="nanopore")
plot_data <- rbind(subset(data,Asm=="nanopore"),newRow1,newRow2)

plot_data$End=as.numeric(levels(plot_data$End))[plot_data$End]

ggplot(aes(y=End, x=Scaffold), data = plot_data) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 17500000, by=500000)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) +
  geom_segment(aes(x=1.8,y=0,xend=1.8,yend=17419488),size=5, color="gray") +
  geom_segment(aes(x=2,y=0,xend=2,yend=246583),size=5, color="gray") +
geom_segment(aes(x=2,y=246583,xend=2,yend=628766),size=5, color="orange") +
geom_segment(aes(x=2,y=628766,xend=2,yend=802382),size=5, color="gray") +
geom_segment(aes(x=2,y=880235,xend=2,yend=955616),size=5, color="orange") +
geom_segment(aes(x=2,y=955653,xend=2,yend=974422),size=5, color="gray") +
geom_segment(aes(x=2,y=1017448,xend=2,yend=1196985),size=5, color="orange") +
geom_segment(aes(x=2,y=1196985,xend=2,yend=1395636),size=5, color="gray") +
geom_segment(aes(x=2,y=1451199,xend=2,yend=1477374),size=5, color="orange") +
geom_segment(aes(x=2,y=1477374,xend=2,yend=1556510),size=5, color="gray") +
geom_segment(aes(x=2,y=1556510,xend=2,yend=1790570),size=5, color="orange") +
geom_segment(aes(x=2,y=1790570,xend=2,yend=2144961),size=5, color="gray") +
geom_segment(aes(x=2,y=2328591,xend=2,yend=2609632),size=5, color="orange") +
geom_segment(aes(x=2,y=2732143,xend=2,yend=3017709),size=5, color="gray") +
geom_segment(aes(x=2,y=3017709,xend=2,yend=3163752),size=5, color="orange") +
geom_segment(aes(x=2,y=3163752,xend=2,yend=3294552),size=5, color="gray") +
geom_segment(aes(x=2,y=3369766,xend=2,yend=3712214),size=5, color="orange") +
geom_segment(aes(x=2,y=3712214,xend=2,yend=3960347),size=5, color="gray") +
geom_segment(aes(x=2,y=4048803,xend=2,yend=4133629),size=5, color="orange") +
geom_segment(aes(x=2,y=4133629,xend=2,yend=4401074),size=5, color="gray") +
geom_segment(aes(x=2,y=4401074,xend=2,yend=5198930),size=5, color="orange") +
geom_segment(aes(x=2,y=5537450,xend=2,yend=5661352),size=5, color="gray") +
geom_segment(aes(x=2,y=5748303,xend=2,yend=6046493),size=5, color="orange") +
geom_segment(aes(x=2,y=6090627,xend=2,yend=6303249),size=5, color="gray") +
geom_segment(aes(x=2,y=6303249,xend=2,yend=6343609),size=5, color="orange") +
geom_segment(aes(x=2,y=6343609,xend=2,yend=6425269),size=5, color="gray") +
geom_segment(aes(x=2,y=6451304,xend=2,yend=6728423),size=5, color="orange") +
geom_segment(aes(x=2,y=6911042,xend=2,yend=8205086),size=5, color="gray") +
geom_segment(aes(x=2,y=8205086,xend=2,yend=8296658),size=5, color="orange") +
geom_segment(aes(x=2,y=8956981,xend=2,yend=9398547),size=5, color="gray") +
geom_segment(aes(x=2,y=9398547,xend=2,yend=9526148),size=5, color="orange") +
geom_segment(aes(x=2,y=9526148,xend=2,yend=9542082),size=5, color="gray") +
geom_segment(aes(x=2,y=9672581,xend=2,yend=10311238),size=5, color="orange") +
geom_segment(aes(x=2,y=10311238,xend=2,yend=11890523),size=5, color="gray") +
geom_segment(aes(x=2,y=11890523,xend=2,yend=11984962),size=5, color="orange") +
geom_segment(aes(x=2,y=12182179,xend=2,yend=12917252),size=5, color="gray") +
geom_segment(aes(x=2,y=12917252,xend=2,yend=13049689),size=5, color="orange") +
geom_segment(aes(x=2,y=14401209,xend=2,yend=16083120),size=5, color="gray") +
geom_segment(aes(x=2,y=16083120,xend=2,yend=17401219),size=5, color="orange") +
  geom_segment(aes(x=2,y=802382,xend=2,yend=880235)) +
  geom_segment(aes(x=2,y=955616,xend=2,yend=955653)) +
  geom_segment(aes(x=2,y=974422,xend=2,yend=1017448)) +
  geom_segment(aes(x=2,y=1395636,xend=2,yend=1451199)) +
  geom_segment(aes(x=2,y=2144961,xend=2,yend=2328591)) +
  geom_segment(aes(x=2,y=2609632,xend=2,yend=2732143)) +
  geom_segment(aes(x=2,y=3294552,xend=2,yend=3369766)) +
  geom_segment(aes(x=2,y=3960347,xend=2,yend=4048803)) +
  geom_segment(aes(x=2,y=5198930,xend=2,yend=5537450)) +
  geom_segment(aes(x=2,y=5661352,xend=2,yend=5748303)) +
  geom_segment(aes(x=2,y=6046493,xend=2,yend=6090627)) +
  geom_segment(aes(x=2,y=6425269,xend=2,yend=6451304)) +
  geom_segment(aes(x=2,y=6728423,xend=2,yend=6911042)) +
  geom_segment(aes(x=2,y=8296658,xend=2,yend=8956981)) +
  geom_segment(aes(x=2,y=9542082,xend=2,yend=9672581)) +
  geom_segment(aes(x=2,y=11984962,xend=2,yend=12182179)) +
  geom_segment(aes(x=2,y=13049689,xend=2,yend=14401209)) +
  geom_segment(aes(x=2,y=17401219,xend=2,yend=17419488)) +
  geom_segment(aes(x=2.2,y=207827,xend=2.2,yend=218997),size=5, color="black") +
  geom_segment(aes(x=2.2,y=6190485,xend=2.2,yend=6278576),size=5, color="black")
 

ref <- data[ which(data$Asm=='ref' & data$Scaffold!="-"), ]
for (row in 1:nrow(ref)) {
  y1 = ref[row,"Start"]
  y2 = ref[row,"End"]
  #line <- paste("geom_segment(aes(x=2,y=",y1,",xend=2,yend=,",y2), sep="")
  plotdata <- plotdata + geom_segment(aes(x=2,y=y1,xend=2,yend=y2))
  cat("geom_segment(aes(x=2,y=",y1,",xend=2,yend=",y2,"))\n",sep="")
}

ref <- data[ which(data$Asm=='ref' & data$Scaffold=="-"), ]
for (row in 1:nrow(ref)) {
  y1 = ref[row,"Start"]
  y2 = ref[row,"End"]
  cat("geom_segment(aes(x=2,y=",y1,",xend=2,yend=",y2,")) +\n",sep="")
}



