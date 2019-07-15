##########################################################################################
#Analysis to examine the relationships among different drivers of biodiversity change
###########################################################################################

#Analysis decisions

#choose ratser data frame file - according to resolution and realm
files <- "output_NP_100_T_RPFirst.RData"#100 km grid terrestrial data
files <- "output_NP_100_M_RPFirst.RData"#100 km grid marine data

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

#get world shapefile
world<-readShapePoly("~/Dropbox/World/world_dissolved.shp")
proj4string(world)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
world<-spTransform(world,CRS(newproj))

####################################################################################

library(reshape2)
library(plyr)

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
#plot chord diagram
#Fig. 1#
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

#Fig. 2 - drawing just climate change correlations
corrMatrix<-cor(mydata@data,method="spearman")
corrMatrixm<-melt(corrMatrix)
corrMatrixm<-subset(corrMatrixm,(Var1%in%c("Temp_trend","VOCC","SST_trend","VOCC_SST")))
corrMatrixm<-subset(corrMatrixm,!(Var1%in%CCvars&Var2%in%CCvars))
corrMatrixm$Var2<-factor(corrMatrixm$Var2,levels=myorder[!myorder%in%CCvars])
corrMatrixm$Var1<-factor(corrMatrixm$Var1)
corrMatrixm$Var1<-factor(corrMatrixm$Var1,levels=levels(corrMatrixm$Var1)[2:1])

ggplot(corrMatrixm)+
  geom_bar(aes(x=Var2,y=value,fill=Var2,alpha=Var1),stat="identity",position="dodge")+
  coord_flip()+
  scale_fill_manual(values=mycols[!myorder%in%CCvars])+
  scale_alpha_manual(values=c(0.8,1))+
  theme_bw()+
  ylim(-0.35,0.35)+
  ylab("Correlation coefficient")+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(axis.text.y = element_text(size=rel(1.5)),
        axis.text.x = element_text(size=rel(1.5)),
        #axis.text.x = element_blank(),axis.ticks.x= element_blank(),
        panel.grid=element_blank(),
        axis.title.x= element_text(size=rel(1.5)),
        axis.title.y= element_blank(),
        legend.position="none")

######################################################################################

#Examining spatial autocorrelation (SOM)

getSpatialAutocorr<-function(df,Var,Increment=10){
  require(ncf)  
  df<-subset(df,!is.na(df[,Var]))
  index = sample(1:length(df$x),2000,replace = F)
  co = correlog(df$x[index],df$y[index],df[,Var][index],increment = Increment,resamp = 1000, 
                latlon = F,na.rm=T,quiet=T)
  co_df<-data.frame(Distance=co[["mean.of.class"]],Corr=co[["correlation"]],Type=Var,P=co[["p"]])
  return(co_df)
}

#run function on the dataset
df<-data.frame(mydataEE@data,mydataEE@coords)
allco<-lapply(names(df)[which(!names(df)%in%c("x","y"))],function(x)getSpatialAutocorr(df,x))
allco<-do.call(rbind,allco)
#split by driver and colour by factor
allco$Driver<-factor(allco$Type,levels=myorder)
allco$Class<-as.character(sapply(allco$Type,getDriverClass))
allco$Class<-factor(allco$Class,levels=c("Climate_change","Human_use","Population",
                                         "Pollution","Alien potential"))
allco<-allco[order(allco$Distance),]

#subset dataset until consistent non-significance for each driver
allco2<-ddply(allco,.(Driver,Class),function(x){
  #first non-sig p
  require(zoo)
  x$rollP<-rollmedian(x$P,5,fill=NA)
  last<-min(which(x$rollP>=0.05&!is.na(x$rollP)))
  x[1:last,]
})

#plotting
png(file = paste0("correlogram_Sig",realm,transformation,".png"),width = 225, height = 550, units = "px")
qplot(Distance,Corr,data=allco,geom="blank",colour=Driver,xlab="Distance in km")+
  theme_bw()+
  stat_smooth(data=allco2,method="loess",alpha=0.9,aes(color=Driver,fill=Driver),se=FALSE)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  facet_wrap(~Class,ncol=1)+
  xlim(0,4000000)+ylim(-0.1,1)+
  ylab("Correlation")+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size=rel(1)))
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
alldataM$Weight<-round(alldataM$weight*100)#get percent cover

#Remove biomes of less interest or not entirely covered by the dataset
alldataM<-subset(alldataM,!BIOME%in%c("Lakes","Rock and Ice","ARCTIC OCEAN","Southern Ocean"))
alldataM$BIOME<-factor(alldataM$BIOME)
alldataM<-subset(alldataM,!is.na(BIOME))
alldataM$Driver<-factor(alldataM$variable,levels=myorder)

#label each variable by its group
alldataM$DriverGroup[alldataM$Driver%in%CCvars]<-"Climate_change"
alldataM$DriverGroup[alldataM$Driver%in%Alien_potential]<-"Alien potential"
alldataM$DriverGroup[alldataM$Driver%in%Pollution]<-"Pollution"
alldataM$DriverGroup[alldataM$Driver%in%HumanUse]<-"Human_use"
alldataM$DriverGroup[alldataM$Driver=="Population"]<-"Human_population"

#get average per cell
alldataM<-ddply(alldataM,.(BIOME,DriverGroup,cell),summarise,
                value = mean(value))

#biome means
avPressureB<-ddply(alldataM,.(DriverGroup,BIOME),summarise,valueB=median(value))
alldataM <- merge(alldataM, avPressureB,by=c("DriverGroup","BIOME"))

#global means
avPressure<-ddply(avPressureB,.(DriverGroup),summarise,value=median(valueB))
alldataM$avPressure<-avPressure$value[match(alldataM$DriverGroup,avPressure$DriverGroup)]

#difference biome means from global ones
alldataM$direction <- as.character(alldataM$valueB>alldataM$avPressure)
alldataM$value<-alldataM$value - alldataM$avPressure

#driver labels sorting
alldataM$DriverGroup<-factor(alldataM$DriverGroup,levels=driverOrder)
if(realm=="T"){
  levels(alldataM$BIOME)<-biomeShort
}
if(realm=="M"){
  alldataM$BIOME<-as.factor(sapply(tolower(as.character(alldataM$BIOME)),simpleCap))
}

#order the biomes by sum of values (total exposure)
sequence<-ddply(alldataM,.(BIOME),summarise,Sum=sum(value))
sequence$BIOME[order(sequence$Sum)]
alldataM$BIOME<-factor(alldataM$BIOME,levels=sequence$BIOME[order(sequence$Sum)])

#plotting
ggplot(alldataM,aes(x=factor(BIOME),y=value))+
  geom_violin(aes(colour=DriverGroup,fill=DriverGroup,alpha=direction))+
  facet_grid(~DriverGroup)+
  scale_alpha_discrete(range=c(0.25,1))+
  scale_fill_manual(values=driverCols)+
  scale_colour_manual(values=driverCols)+
  coord_flip()+theme_bw()+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(breaks=c(-0.6,-0.3,0,0.3,0.6))+
  facet_grid(~DriverGroup)+
  ylab("Deviation from global average")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=rel(1.5),angle=90),
        axis.text.y = element_text(size=rel(1.5)),
        axis.text.x = element_text(size=rel(1.2)),
        axis.title.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        #axis.line = element_line(colour = "black"))
        legend.position="none")

#######################################################################################  

###################
#Cluster analysis#
##################

clustdata <- mydata@data

#apply cluster analysis
library(fpc)
library(cluster)
out<-lapply(2:10,function(x){
  pam(clustdata, x, metric="euclidean")
})

#check metrics by the number of clusters:
#dissimilarity
#average dissimilarity between the observations in the cluster and the cluster's medoid
q1<- qplot(2:10,sapply(out,function(x)median(x$clusinfo[,2])))+ggtitle("dissim")
#separation
#minimal dissimilarity between an observation of the cluster and an observation of another cluster
q2 <- qplot(2:10,sapply(out,function(x)median(x$clusinfo[,5])))+ggtitle("sep")
#silhouette widths
q3 <- qplot(2:10,sapply(out,function(x)median(x$silinfo$clus.avg.width)))+ggtitle("med")

#plot all
library(cowplot)
plot_grid(q1,q2,q3)


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

##########################################################################

#Plotting global maps

#assign 1 to cells in the top 10% of non-zero vaues
mydataP<-mydataR
mydataP<-spTransform(mydataP,projection(refProj))

mydataP@data<-data.frame(sapply(names(mydataP@data),function(x){
  var<-mydataP@data[,x]
  mydataP@data[,x]<-ifelse(mydataP@data[,x]<=0,0,mydataP@data[,x])
  quant<-quantile(mydataP@data[,x][mydataP@data[,x]>0],0.9)
  mydataP@data[,x]<-ifelse(mydataP@data[,x]<=quant,0,1)
}))

colSums(mydataP@data)
nrow(mydata@data)

#function to sum across rasters within a driver group
myfun<-function(mydataP){
  
  #get number of rasters within driver  
  N=length(names(mydataP))
  
  #find hotspot areas 
  rstack<-stack()
  for(i in 1:N){
    out<-rasterize(mydataP,refProj,field=names(mydataP)[i])
    rstack<-stack(list(rstack,out))
  }
  out<-calc(rstack,fun=sum)
  return(out)
}

#create dummy variable for later use
mydataP@data$dummy<-rep(0,nrow(mydataP))

#code for terrestrial data sets
outT_TC<-myfun(mydataP[,CCvars])
#for human use
outT_HU<-myfun(mydataP[,c("Crop_trend","Cropland","Forest_loss","Pasture_trend","Urban","Urban_trend")])
#for human population
outT_Pop<-myfun(mydataP[,c("dummy","Population")])
#for pollution
outT_P<-myfun(mydataP[,c("Fertilizer_app","Pesticides","N_deposition")])
#for Alien potential
outT_I<-myfun(mydataP[,c("dummy","Accessibility")])
outSum<-stack(outT_TC,outT_P,outT_HU,outT_Pop,outT_I)
save(outSum,file="outSum.RData")

#code for marine datasets
outM_TC<-myfun(mydataP[,c(CCvars)])
#for Alien potential
outM_I<-myfun(mydataP[,c("dummy","Port_volume")])
#for human population
outM_Pop<-myfun(mydataP[,c("dummy","Population")])
#for pollution
outM_P<-myfun(mydataP[,c("Fertilizer","Ocean_poll","Inorganic" )])
#for human use
outM_HU<-myfun(mydataP[,c("Artisanal_fish","Demersalfish_Destr","Demersalfish_HighBycatch",
                          "Demersalfish_LowBycatch","Pelagicfish_HighBycatch",
                          "Pelagicfish_LowBycatch")])
outSum_M<-stack(outM_TC,outM_P,outM_HU,outM_Pop,outM_I)
save(outSum,file="outSum_M.RData")

#Plotting combined map
load("outSum.RData")
load("outSum_M.RData")

library(RColorBrewer)
library(scales)
show_col(brewer.pal(9,"Oranges"))

#plotting 5 in a row
png("global_maps.png", 550, 500)
par(mfrow=c(3,2))

#merge terrestrial and marine maps for each driver group
outB_TC<-merge(outSum[[1]],outSum_M[[1]])
plot(outB_TC,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[c(3,5,6,9)]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.1)

outB_HU<-merge(outSum_M[[3]],outSum[[3]])
plot(outB_HU,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[c(3,4,6,7,8,9)]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.1)

outB_Pop<-merge(outSum[[4]],outSum_M[[4]])
plot(outB_Pop,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[9]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.1)

outB_P<-merge(outSum[[2]],outSum_M[[2]])
plot(outB_P,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[c(3,6,9)]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.1)

outB_I<-merge(outSum[[5]],outSum_M[[5]])
plot(outB_I,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[9]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.1)

#summed map
outB_sum<-stack(outB_TC,outB_HU,outB_Pop,outB_P,outB_I)
outB_sum<-calc(outB_sum,fun=sum)
plot(outB_sum,col=c(col2hex("gray99"),brewer.pal(9,"Oranges")[c(2,3,4,5,7,8,9)]),axes=FALSE,legend=TRUE, asp = 1)
plot(world,add=T,border=col2hex("gray20"),lwd=0.5)
dev.off()

###########################################################################
