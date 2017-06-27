##########################################################################################
#Analysis to examine the relationships among different drivers of biodiversity change
###########################################################################################

#Analysis decisions

#choose ratser data frame file - according to resolution and projection
files<-"output_all_Eckman_100km.RData"

#decide whether to do analysis for terrestrial or for marine
realm<-"T"
realm<-"M"

#choose transformation type of the raster values
transformation="rank0"
transformation="log"

###############################################################################

#Producing reference raster grid

library(raster)

#set equal area projection 
newproj<-"+proj=eck4 +datum=WGS84"
newres=100000

#create 1 degree grid and convert it into Eckert IV 100 km reference grid
ref<-extent(-180, 180, -90, 90)
ref<-raster(ref)
res(ref)<-1
values(ref)<-1
projection(ref)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
refProj<-projectRaster(ref, crs=newproj, res=newres,over=T)
plot(refProj)

##############################################################################

#run global data procesing script
source('~/Desktop/processing.R', echo=TRUE)

########################################################################################### 

#get world dissolved
world<-readShapePoly("~/Dropbox/World/world_dissolved.shp")
proj4string(world)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
world<-spTransform(world,CRS(newproj))


####################################################################################

######################
#Correlation analysis#
######################

#get correlation matrix
corrMatrix<-cor(mydata@data,method="spearman")

library(circlize)
corrMatrix[upper.tri(corrMatrix)] <- NA
#set everything less than 0.7 as 0 (for transparency, see later)
corrMatrix[corrMatrix<0.7] <- 0.0

#melt matrix and remove identity correlations
corrMatrixm<-melt(corrMatrix)
corrMatrixm<-subset(corrMatrixm,!is.na(value))
corrMatrixm<-subset(corrMatrixm,value!=1)

#specific colour of strong correlation links
corrMatrixm$Colour[corrMatrixm$value!=0.1]<-col2hex("grey70")
#shade out weak links
corrMatrixm$Colour[corrMatrixm$value==0.0]<-"#FFFFFF00"

#####################
#plot chord diagram#
####################
png(file = paste0("chordplot_",realm,transformation,"0.7.png"),width = 500, height = 500, units = "px")

chordDiagram(corrMatrixm,symmetric = FALSE,
             transparency=0.5,
             col=corrMatrixm$Colour,
             grid.col=rev(mycols),
             order=rev(myorder),
             annotationTrack = "grid", preAllocateTracks = 1)
#change label direction
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, cex=0.6,facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.2, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

######################################################################################

#modified t.test to check significance of pair-wise correlations after accounting 
#for spatial autocorrelation

library(SpatialPack)
out<-apply(corrMatrixm,1,function(x){
  modified.ttest(x=mydataEE@data[,x["Var1"]],y=mydataEE@data[,x["Var2"]],coords=mydataEE@coords)})

#Obtaining the corrected p-values
correctedCors<-data.frame(Var1=corrMatrixm$Var1,Var2=corrMatrixm$Var2,value=corrMatrixm$value,P=sapply(out,function(x)x$p.value))

#terrestrial
subset(correctedCors,value>=0.7)#all significant
#marine
subset(correctedCors,value>=0.7)# all significant

#######################################################################################

###################
#Biome differences#
###################

#get data frame showing biome overlap for each grid cell
biomeCov<-ddply(biomeCov,.(cell,BIOME,x,y),summarise,weight=sum(weight))
alldata<-data.frame(mydataEE@data,mydataEE@coords)
alldata<-merge(biomeCov,alldata,by=c("x","y"))
alldataM<-melt(alldata,id=c("x","y","BIOME","cell","weight"))

#Remove biomes of less interest or not entirely covered by the dataset
alldataM<-subset(alldataM,!BIOME%in%c("Lakes","Rock and Ice","ARCTIC OCEAN","Southern Ocean"))
alldataM$BIOME<-factor(alldataM$BIOME)
alldataM<-subset(alldataM,!is.na(BIOME))
alldataM$Driver<-factor(alldataM$variable,levels=myorder)

#For each biome and driver, get weighted average rank
alldataMa<-ddply(alldataM,.(BIOME,Driver),summarise,
                 sd=sd(rep(value,times=Weight)),
                 value=mean(rep(value,times=Weight)))


#Difference by the global average
avPressure<-ddply(alldataMa,.(Driver),summarise,value=median(value))
alldataMa$avPressure<-avPressure$value[match(alldataMa$Driver,avPressure$Driver)]
alldataMa$value<-(alldataMa$value-alldataMa$avPressure)
alldataMa<-subset(alldataMa,value>0)#just plot the positive deviations

#Cleaning for presentation
alldataMa<-alldataMa[order(as.numeric(alldataMa$Driver)),]
alldataMa$BIOME<-as.factor(sapply(tolower(as.character(alldataMa$BIOME)),simpleCap))

#terrestrial labels
levels(alldataMa$BIOME)<-biomeShort

#order the biome variable by the number of climate change impacts it has, and then sum
sequence<-ddply(alldataMa,.(BIOME),summarise,Sum=sum(value))
alldataMa$BIOME<-factor(alldataMa$BIOME,levels=sequence$BIOME[order(sequence$Sum)])

#Setting order of panels
alldataMa$DriverF<-as.factor(alldataMa$Driver)

#shorten the fishing ones
levels(alldataMa$DriverF)[6:10]<-c("Pelagic_LowBycatch","Pelagic_HighBycatch",
                                   "Demersal_Destr","Demersal_LowBycatch",
                                   "Demersal_HighBycatch")

alldataMa$DriverF<-factor(alldataMa$DriverF,levels=rev(levels(alldataMa$DriverF)))

#dot map plot
png(file = paste0("sepedchart2square",realm,transformation,".png"),width = 1200, height = 850, units = "px")
ggplot(alldataMa,aes(x=BIOME,y=value))+
  geom_bar(aes(fill=Driver),stat="identity")+
  scale_fill_manual(values=mycols)+
  geom_errorbar(aes(x=BIOME,ymax=value+sd,ymin=value,colour=Driver))+
  scale_colour_manual(values=mycols)+
  coord_flip()+theme_bw()+
  facet_grid(~DriverF)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=rel(1.5),angle=90),
        axis.text.y = element_text(size=rel(1.5)),
        axis.text.x = element_blank(),axis.ticks.x= element_blank(),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        legend.position="none")
dev.off()

#######################################################################################  

###################
#Cluster analysis#
##################

clustdata <- mydata@data

#######################################  
#first using PCA to condense the data##
#######################################

head(clustdata)
scores<-data.frame(Population=clustdata[,Population],Invasion=clustdata[,Invasions])

#need to run the following code line by line
for(var in c(CCvars,HumanUse,Pollution)){
  fit <- princomp(clustdata[,var], cor=TRUE)
  print(loadings(fit,cutoff=0.0))
  summary(fit)
  temp<-data.frame(fit$scores[,1:2])
  #names(temp)<-c(paste0(var[1],"1"),paste0(var[1],"2"))
  names(temp)<-c(paste0(var[1],"1"),paste0(var[1],"2"),paste0(var[1],"3"))
  scores<-data.frame(scores,temp)
}

#apply cluster analysis
library(fpc)
library(cluster)
out<-lapply(2:10,function(x){
  pam(scores, x, metric="manhattan")
})

save(out,file=paste0("out_PCAcluster",realm,transformation,"2.RData"))

#check metrics by the number of clusters

#dissimilarity
qplot(2:10,sapply(out,function(x)mean(x$clusinfo[,2])))

#silhouette width
plot(2:10,sapply(out,function(x)median(x$silinfo$clus.avg.width)))
plot(2:10,sapply(out,function(x)mean(x$silinfo$clus.avg.width)))
plot(2:10,sapply(out,function(x)min(x$silinfo$clus.avg.width)))

#Selecting number of clusters

#for terrestrial
clusterNumber<-5 

#for marine
clusterNumber<-6

#Negative widths reassigned to neighbours
outSI<-data.frame(out[[clusterNumber-1]]$silinfo$widths)
outSI$rN<-as.numeric(row.names(outSI))
outSI<-outSI[order(outSI$rN),]
outSI$cluster[outSI$sil_width<0]<-outSI$neighbor[outSI$sil_width<0]

#final classification
fit<-out[[clusterNumber-1]]
grp<-outSI$cluster

#####################
#Plotting the legend#
#####################

#to work out the meaning of each cluster group
clustdataDF<-data.frame(clustdata,grp)
clustDFmelted<-melt(clustdataDF,id="grp")

#what is the highest pressure in each group?
centre<-ddply(clustdataDF,.(grp),function(x)colMeans(x))[,1:16]
selectVars<-apply(centre,1,function(x)names(centre)[x==max(x)])
colsT2<-mycols[match(selectVars,myorder)]

#Plotting the legend
clustDFmelted$variable<-factor(clustDFmelted$variable,levels=myorder)
clustDFmelted$grp<-as.factor(clustDFmelted$grp)

png(file = paste0("clustermap_BarChartSide",realm,transformation,".png"),width = 500, height = 275, units = "px")
ggplot(data=clustDFmelted)+
  geom_boxplot(aes(x=variable,y=value,colour=variable,fill=variable),
               outlier.shape=NA,outlier.size = 0, coef = 0)+
  facet_wrap(~grp,ncol=16)+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        legend.position="none",
        axis.text.y = element_text(size=14)) 
dev.off()


#########################
#Plotting the cluster map#
#########################

alldata<-data.frame(mydataEE@data,mydataEE@coords)
alldata$groups<-grp
coordinates(alldata)<-c("x","y")
proj4string(alldata)<-newproj

#create a raster from the data
alldata$groups<-as.numeric(alldata$groups)
r<-rasterize(alldata,refProj,field="groups")

#Slighly smooth the groups spatially by taking the mode cluster group per 3 x 3 neighbour
w=matrix(1,nrow=3,ncol=3)
r <- focal(r, w=w,fun=modal,na.rm=T)

#create smoother continental edges for presentation
rd <- disaggregate(r, fact=c(4, 4))

#Masking 
library(maptools)
data(wrld_simpl)
wrld_simpl<-spTransform(wrld_simpl,CRS(projection(refProj)))
if (realm=="T"){
  rd<-mask(rd,wrld_simpl,inverse=F)
}else if (realm=="M"){
  rd<-mask(rd,wrld_simpl,inverse=T)  
}

#plot raster in ggplot
g<-gplot(rd) + geom_tile(aes(fill = factor(value)))+
  scale_fill_manual(values=colsT2,na.value="white") +
  coord_equal()+theme(legend.position="none")

#extract the realm boundaries
worldP<-world
worldP@data$id = rownames(worldP@data)
worldP = fortify(worldP, region="id")

#adding the world border
pp <- rasterToPolygons(refProj, dissolve=TRUE)
outline <- fortify(pp)
png(file = paste0("clustermap_border_",realm,transformation,clusterNumber,".png"),width = 900, height = 700, units = "px")
g+geom_polygon(data=worldP,aes(x=long,y=lat,group=group),fill="NA",colour="black",size=0.25)+
  geom_path(aes(x = long, y = lat), data = outline, size=0.25, colour="black")
dev.off()  

###########################################################################

