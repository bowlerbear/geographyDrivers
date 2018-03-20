###########################################################################################
#Processsing script to covert the raster data frame into a form that is ready for analysis
###########################################################################################

#get libraries well need
library(raster)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(rasterVis)
library(ncdf4)
library(plyr)
library(maptools)
library(maps)
library(reshape2)

######################################################################

#function to read each file
rawData<-lapply(files,function(x){
  myfile<-load(x)
  myfile<-get(myfile)
  myfile<-do.call(rbind,myfile)
  myfile$File<-x
  return(myfile)
})

#combine all raster data frames
rawData<-do.call(rbind,rawData)

#########################################################################################

#add realm info to our data frame
setwd("~/Documents/sChange/driver_maps")
myrasters<-read.csv("myrasters.csv",sep=";",as.is=T)
rawData$Realm<-myrasters$Realm[match(rawData$Type,myrasters$name)]

###########################################################################################

#select data frames to be used

#terrestrial
if(realm=="T"){
  rawData<-subset(rawData,Realm=="T")
  
}else if(realm=="M"){
  rawData<-subset(rawData,Realm=="M")
}

###############################################################################

#plot where we have data

mydata<-rawData
mydata<-subset(mydata,!is.na(mydata$Value))
coordinates(mydata)<-c("x","y")
proj4string(mydata)<-CRS(projection(refProj))
mydata<-spTransform(mydata,"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
mydata@data$x<-mydata@coords[,1]
mydata@data$y<-mydata@coords[,2]

ggplot(mydata@data)+geom_point(aes(x=x,y=y),size=rel(0.5))+facet_wrap(~Type)

############################################################################

#how to cope with missing data?

#assumed that it is 0 - that seems reasonable looking at where the missing values are

#terrestrial
if(realm=="T"){
  #most missing data in very north of northern hemisphere assume as zero
  rawData$Value[rawData$Type=="Pesticides"&is.na(rawData$Value)]<-0.00
  rawData$Value[rawData$Type=="Livestock"&is.na(rawData$Value)]<-0.00
  rawData$Value[rawData$Type=="Cropland"&is.na(rawData$Value)]<-0.00
  rawData$Value[rawData$Type=="Pasture"&is.na(rawData$Value)]<-0.00
  rawData$Value[rawData$Type=="Pop_dens"&is.na(rawData$Value)]<-0.00
  
}else if (realm=="M"){
  rawData$Value[rawData$Type%in%c("Artisanal_fish","ArtisanalFish")&is.na(rawData$Value)]<-0.00
}

############################################################################################

#use just grid cells where we have a complete dataset
#crop to a constant extent

rawData<-dcast(rawData,x+y~Type,value.var="Value")
mydata<-rawData[complete.cases(rawData),]
refCrop<-crop(ref,extent(-179,179,-58,78))
out<-projectExtent(refCrop,crs=newproj)
mydata<-subset(mydata,x>out@extent@xmin&x<out@extent@xmax)
mydata<-subset(mydata,y>out@extent@ymin&y<out@extent@ymax)

#########################################################################################

#where do we have data now?

coordinates(mydata)<-c("x","y")
proj4string(mydata)<-CRS(newproj)
mydataEE<-mydata#equal area projection
mydata<-spTransform(mydata,"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

#plotting
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
plot(mydata,col="red",add=T)

#save as grid
mygrid<-rasterize(mydataEE,refProj,field=names(mydataEE)[1])
r<-calc(mygrid,fun=function(x)ifelse(!is.na(x),1,NA))
writeRaster(r,file=paste0("grid_",realm,newres,".tif"),format="GTiff",overwrite=T)

#######################################################################################

#Plotting the distribution of values of each data layer

library(gridExtra)

p<-lapply(names(mydata@data),function(x){
  qplot(x=mydata@data[,x],data=mydata@data,geom="histogram")+
    ggtitle(x)+
    xlab("Raster units")+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
})


do.call(grid.arrange, p)

#####################################################################################

#Make transformations so that larger values mean more of that pressure

if(realm=="T"){
  mydata@data$Accessibility<-1/mydata@data$Accessibility
  mydata@data$Forest_loss<-mydata@data$Forest_loss*-1
}

#####################################################################################

#Raster value transformation 

if(transformation=="rank0"){
  mydata@data<-data.frame(sapply(names(mydata@data),function(x){
    temp<-rank(unique(mydata@data[,x]))
    mydata@data[,x]<-temp[match(mydata@data[,x],unique(mydata@data[,x]))]
  }))
  #also scale between 0 and 1
  mydata@data<-data.frame(sapply(names(mydata@data),function(x){
    mydata@data[,x]<-(mydata@data[,x]-min(mydata@data[,x]))/(max(mydata@data[,x])-min(mydata@data[,x]))
  }))
}else if (transformation=="log"){
  mydata@data<-data.frame(sapply(names(mydata@data),function(x){
    #bound upper and lower values
    upper<-quantile(mydata@data[,x],0.975)
    mydata@data[,x]<-sapply(mydata@data[,x],function(j)ifelse(j>upper,upper,j))
    lower<-quantile(mydata@data[,x],0.025)
    mydata@data[,x]<-sapply(mydata@data[,x],function(j)ifelse(j<lower,lower,j))
  }))
  mydata@data<-data.frame(sapply(names(mydata@data),function(x){
    if(all(mydata@data[,x]>=0)){
      mydata@data[,x]<-log10(mydata@data[,x]+1)
    }else{
      mydata@data[,x]<-mydata@data[,x]
    }
  }))
  #also scale between 0 and 1
  mydata@data<-data.frame(sapply(names(mydata@data),function(x){
    mydata@data[,x]<-(mydata@data[,x]-min(mydata@data[,x]))/(max(mydata@data[,x])-min(mydata@data[,x]))
  }))
}}


#look at summary of the raster values as a check
summary(mydata@data)
dim(mydata@data)
p<-lapply(names(mydata@data),function(x){
  qplot(x=mydata@data[,x],data=mydata@data,geom="histogram")+ggtitle(x)
})

do.call(grid.arrange, p)

#####################################################################################

#plot data again
mydataP<-mydata
mydataP@data$x<-mydataP@coords[,1]
mydataP@data$y<-mydataP@coords[,2]
mydataPmelt<-melt(mydataP@data,id=c("x","y"))
ggplot(mydataPmelt)+geom_point(aes(x=x,y=y,color=value),size=rel(0.2))+
    facet_wrap(~variable)+
    scale_color_gradient(low="white",high="steelblue")+
    xlab("Longitude (E)")+
    ylab("Latitude (N)")+
    theme(strip.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_rect(colour="white"))

#########################################################################################

#Also make an equal area projected data set
mydataEE<-spTransform(mydata,CRS(newproj))

################################################################################################

#Setting plot order and colours

library(RColorBrewer)
library(gplots)

if(realm=="T"){
  myorder=order=c("Accessibility",
                  "Pesticides", "Fertilizer_app","N_deposition",
                  "Population",
                  "Pasture_trend","Cropland","Crop_trend","Urban","Urban_trend","Forest_loss",
                  "Aridity_change","Temp_divergence","Extreme_trends","VOCC","Temp_change")
  
  biomeAbrevs<-c("BF","D","FS","L","M","MS","MGS","RI",
                 "TBMF","TeCF","TeGS","TrCF","TDBF","TrGS",
                 "TMBF","T")
  
  biomeShort<-c("Boreal forest","Desert","Flooded savanna","Mangrove",
                "Mediterranean scrub","Montane grass/shrubland","Temperate broadleaf/mixed forest",
                "Temperate coniferous forest","Temperate grass/shrubland","(Sub)Tropical coniferous forest",
                "(Sub)Tropical dry broadleaf forest","(Sub)Tropical grass/shrubland","(Sub)Tropical moist broadleaf forest",
                "Tundra")
  
  Climate_change<-CCvars<-c("Temp_change","Extreme_trends","VOCC",
                            "Temp_divergence","Aridity_change")
  Human1<-c("Fertilizer_app","N_deposition","Pesticides","Accessibility","Urban_trend","Cropland")
  Human2<-c("Pasture_trend","Forest_loss","Crop_trend")
  Human_use<-HumanUse<-c("Urban_trend","Urban","Cropland","Pasture_trend","Forest_loss","Crop_trend")
  Invasions<-"Accessibility"
  Pollution<-c("Fertilizer_app","N_deposition","Pesticides")
  
}else if (realm=="M"){
  myorder<-rev(c("SST_change","VOCC_SST","SST_extremes","SST_divergence","Ocean_acid",
                 "Artisanal_fish","Demersalfish_HighBycatch","Demersalfish_LowBycatch", 
                 "Demersalfish_Destr","Pelagicfish_HighBycatch","Pelagicfish_LowBycatch",
                 "Population",
                 "Fertilizer","Inorganic","Ocean_poll",
                 "Port_volume"))
  
  biomeAbrevs<-c("A","CIP","EIP","IO","NAO","NPO",
                 "SAO","SPO",
                 "TAu","TNA","TNP","TSAm","TSAf","TAt",
                 "TEP","WIP")
  
  Climate_change<-CCvars<-c("Ocean_acid","SST_change","SST_extremes","SST_divergence","VOCC_SST")
  Human1<-c("Artisanal_fish","Port_volume","Fertilizer","Inorganic")
  Human_use<-HumanUse<-c("Artisanal_fish","Demersalfish_Destr","Demersalfish_HighBycatch","Demersalfish_LowBycatch",
                         "Pelagicfish_HighBycatch","Pelagicfish_LowBycatch")
  Invasions<-"Port_volume"
  Pollution<-c("Fertilizer","Inorganic","Ocean_poll") 
}

mycols<-c(mycols<-col2hex("lightblue4"),
          brewer.pal(9,"Greys")[5:7],
          brewer.pal(9,"Purples")[5],
          brewer.pal(9,"Blues")[2:7],
          brewer.pal(9,"OrRd")[3:7])

driverOrder<-c("Climate_change","Human_use","Population","Pollution","Invasions")
driverCols<-c(brewer.pal(9,"OrRd")[7],brewer.pal(9,"Blues")[7],brewer.pal(9,"Purples")[5],
              brewer.pal(9,"Greys")[7],brewer.pal(11,"BrBG")[5])


####################################################################################################


