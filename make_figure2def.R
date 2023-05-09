library("phyloseq")
library("HybridMicrobiomes")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo17 Leaf 16S Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data<-read.csv(paste0('Maize/genus_bacteria.csv'), row.names=1,header= TRUE)
plantdata<-read.table('Maize/planting_design2017.txt',header=TRUE)
microbedata<-read.table('Maize/DNA_HYB.txt',header=TRUE)

pureB73<-plantdata$PlantID[intersect(which(plantdata$F_parent=='B73'), which(plantdata$M_parent=='B73'))]
pureMo17 <-plantdata$PlantID[intersect(which(plantdata$F_parent=='Mo17'), which(plantdata$M_parent=='Mo17'))]
hybridB73Mo17<-plantdata$PlantID[intersect(which(plantdata$F_parent=='B73'), which(plantdata$M_parent=='Mo17'))]

pureB73microbe<-which(microbedata$PlantID %in% pureB73)
pureMo17microbe<-which(microbedata$PlantID %in% pureMo17)
hybridB73Mo17microbe<-which(microbedata$PlantID %in% hybridB73Mo17)

plant_part<-which(microbedata$Type=='leaf')

pureB73leaf<-microbedata$SampleID[intersect(pureB73microbe,plant_part)]
pureMo17leaf<-microbedata$SampleID[intersect(pureMo17microbe,plant_part)]
hybridB73Mo17leaf<-microbedata$SampleID[intersect(hybridB73Mo17microbe,plant_part)]

data<-t(data)
matrixpureB73leaf<-which(rownames(data) %in% pureB73leaf)
matrixpureMo17leaf<-which(rownames(data) %in% pureMo17leaf)
matrixhybridB73Mo17leaf<-which(rownames(data) %in% hybridB73Mo17leaf)

OTUB73leaf<-data[matrixpureB73leaf,]
OTUMo17leaf<-data[matrixpureMo17leaf,]
OTUB73Mo17leaf<-data[matrixhybridB73Mo17leaf,]

#Remove all samples with fewer than this number of reads (note that this number should be much smaller if you are removing environmental microbes)
#Based on samples, I suggest 10000 when not removing controls and 1500 when removing controls
readthreshold<-5000

datarfullB73Mo17leaf<-t(rbind(OTUB73leaf,OTUMo17leaf,OTUB73Mo17leaf))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo17leaf)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo17leaf<-datarfullB73Mo17leaf[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo17leaf<-datarfullB73Mo17leaf
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo17leaf = otu_table(datarB73Mo17leaf, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo17leaf[,1])){
  temp<-strsplit(rownames(datarB73Mo17leaf)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo17leaf)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo17leaf<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo17leaf_physeq = phyloseq(OTUB73Mo17leaf, TAXB73Mo17leaf)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo17leaf_physeq)


#Get species/control-type information
B73Mo17leafno<-c()
for (k in 1:length(sample_name)){
if (sample_name[k] %in% pureB73leaf){
  B73Mo17leafno<-c(B73Mo17leafno,1)
}
  else if (sample_name[k] %in% pureMo17leaf){
    B73Mo17leafno<-c(B73Mo17leafno,3)
  }
  else if (sample_name[k] %in% hybridB73Mo17leaf){
    B73Mo17leafno<-c(B73Mo17leafno,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo17leaf_bs<-FourHbootstrap(B73Mo17leaf_physeq,B73Mo17leafno,0.5,500,10)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo17 Rhizosphere 16S Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plant_part<-which(microbedata$Type=='rhizosphere')

pureB73rhizosphere<-microbedata$SampleID[intersect(pureB73microbe,plant_part)]
pureMo17rhizosphere<-microbedata$SampleID[intersect(pureMo17microbe,plant_part)]
hybridB73Mo17rhizosphere<-microbedata$SampleID[intersect(hybridB73Mo17microbe,plant_part)]

matrixpureB73rhizosphere<-which(rownames(data) %in% pureB73rhizosphere)
matrixpureMo17rhizosphere<-which(rownames(data) %in% pureMo17rhizosphere)
matrixhybridB73Mo17rhizosphere<-which(rownames(data) %in% hybridB73Mo17rhizosphere)

OTUB73rhizosphere<-data[matrixpureB73rhizosphere,]
OTUMo17rhizosphere<-data[matrixpureMo17rhizosphere,]
OTUB73Mo17rhizosphere<-data[matrixhybridB73Mo17rhizosphere,]

datarfullB73Mo17rhizosphere<-t(rbind(OTUB73rhizosphere,OTUMo17rhizosphere,OTUB73Mo17rhizosphere))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo17rhizosphere)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo17rhizosphere<-datarfullB73Mo17rhizosphere[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo17rhizosphere<-datarfullB73Mo17rhizosphere
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo17rhizosphere = otu_table(datarB73Mo17rhizosphere, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo17rhizosphere[,1])){
  temp<-strsplit(rownames(datarB73Mo17rhizosphere)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo17rhizosphere)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo17rhizosphere<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo17rhizosphere_physeq = phyloseq(OTUB73Mo17rhizosphere, TAXB73Mo17rhizosphere)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo17rhizosphere_physeq)


#Get species/control-type information
B73Mo17rhizosphereno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73rhizosphere){
    B73Mo17rhizosphereno<-c(B73Mo17rhizosphereno,1)
  }
  else if (sample_name[k] %in% pureMo17rhizosphere){
    B73Mo17rhizosphereno<-c(B73Mo17rhizosphereno,3)
  }
  else if (sample_name[k] %in% hybridB73Mo17rhizosphere){
    B73Mo17rhizosphereno<-c(B73Mo17rhizosphereno,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo17rhizosphere_bs<-FourHbootstrap(B73Mo17rhizosphere_physeq,B73Mo17rhizosphereno,0.5,500,10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo17 Leaf ITS Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data<-read.csv(paste0('Maize/fungi.csv'), row.names=1,header= TRUE)

data<-t(data)

matrixpureB73leaf<-which(rownames(data) %in% pureB73leaf)
matrixpureMo17leaf<-which(rownames(data) %in% pureMo17leaf)
matrixhybridB73Mo17leaf<-which(rownames(data) %in% hybridB73Mo17leaf)

OTUB73leaffungi<-data[matrixpureB73leaf,]
OTUMo17leaffungi<-data[matrixpureMo17leaf,]
OTUB73Mo17leaffungi<-data[matrixhybridB73Mo17leaf,]

datarfullB73Mo17leaffungi<-t(rbind(OTUB73leaffungi,OTUMo17leaffungi,OTUB73Mo17leaffungi))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo17leaffungi)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo17leaffungi<-datarfullB73Mo17leaffungi[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo17leaffungi<-datarfullB73Mo17leaffungi
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo17leaffungi = otu_table(datarB73Mo17leaffungi, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo17leaffungi[,1])){
  temp<-strsplit(rownames(datarB73Mo17leaffungi)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo17leaffungi)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo17leaffungi<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo17leaffungi_physeq = phyloseq(OTUB73Mo17leaffungi, TAXB73Mo17leaffungi)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo17leaffungi_physeq)


#Get species/control-type information
B73Mo17leaffungino<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73leaf){
    B73Mo17leaffungino<-c(B73Mo17leaffungino,1)
  }
  else if (sample_name[k] %in% pureMo17leaf){
    B73Mo17leaffungino<-c(B73Mo17leaffungino,3)
  }
  else if (sample_name[k] %in% hybridB73Mo17leaf){
    B73Mo17leaffungino<-c(B73Mo17leaffungino,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo17leaffungi_bs<-FourHbootstrap(B73Mo17leaffungi_physeq,B73Mo17leaffungino,0.5,500,10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo17 Rhizosphere ITS Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matrixpureB73rhizosphere<-which(rownames(data) %in% pureB73rhizosphere)
matrixpureMo17rhizosphere<-which(rownames(data) %in% pureMo17rhizosphere)
matrixhybridB73Mo17rhizosphere<-which(rownames(data) %in% hybridB73Mo17rhizosphere)

OTUB73rhizospherefungi<-data[matrixpureB73rhizosphere,]
OTUMo17rhizospherefungi<-data[matrixpureMo17rhizosphere,]
OTUB73Mo17rhizospherefungi<-data[matrixhybridB73Mo17rhizosphere,]

datarfullB73Mo17rhizospherefungi<-t(rbind(OTUB73rhizospherefungi,OTUMo17rhizospherefungi,OTUB73Mo17rhizospherefungi))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo17rhizospherefungi)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo17rhizospherefungi<-datarfullB73Mo17rhizospherefungi[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo17rhizospherefungi<-datarfullB73Mo17rhizospherefungi
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo17rhizospherefungi = otu_table(datarB73Mo17rhizospherefungi, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo17rhizospherefungi[,1])){
  temp<-strsplit(rownames(datarB73Mo17rhizospherefungi)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo17rhizospherefungi)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo17rhizospherefungi<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo17rhizospherefungi_physeq = phyloseq(OTUB73Mo17rhizospherefungi, TAXB73Mo17rhizospherefungi)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo17rhizospherefungi_physeq)


#Get species/control-type information
B73Mo17rhizospherefungino<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73rhizosphere){
    B73Mo17rhizospherefungino<-c(B73Mo17rhizospherefungino,1)
  }
  else if (sample_name[k] %in% pureMo17rhizosphere){
    B73Mo17rhizospherefungino<-c(B73Mo17rhizospherefungino,3)
  }
  else if (sample_name[k] %in% hybridB73Mo17rhizosphere){
    B73Mo17rhizospherefungino<-c(B73Mo17rhizospherefungino,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo17rhizospherefungi_bs<-FourHbootstrap(B73Mo17rhizospherefungi_physeq,B73Mo17rhizospherefungino,0.5,500,10)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xCML103 Leaf 16S Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data<-read.csv(paste0('Maize/genus_bacteria.csv'), row.names=1,header= TRUE)
plantdata<-read.table('Maize/planting_design2017.txt',header=TRUE)
microbedata<-read.table('Maize/DNA_HYB.txt',header=TRUE)

pureCML103 <-plantdata$PlantID[intersect(which(plantdata$F_parent=='CML103'), which(plantdata$M_parent=='CML103'))]
hybridB73CML103<-plantdata$PlantID[intersect(which(plantdata$F_parent=='B73'), which(plantdata$M_parent=='CML103'))]

pureCML103microbe<-which(microbedata$PlantID %in% pureCML103)
hybridB73CML103microbe<-which(microbedata$PlantID %in% hybridB73CML103)

plant_part<-which(microbedata$Type=='leaf')

pureCML103leaf<-microbedata$SampleID[intersect(pureCML103microbe,plant_part)]
hybridB73CML103leaf<-microbedata$SampleID[intersect(hybridB73CML103microbe,plant_part)]

data<-t(data)
matrixpureB73leaf<-which(rownames(data) %in% pureB73leaf)
matrixpureCML103leaf<-which(rownames(data) %in% pureCML103leaf)
matrixhybridB73CML103leaf<-which(rownames(data) %in% hybridB73CML103leaf)

OTUB73leaf<-data[matrixpureB73leaf,]
OTUCML103leaf<-data[matrixpureCML103leaf,]
OTUB73CML103leaf<-data[matrixhybridB73CML103leaf,]

datarfullB73CML103leaf<-t(rbind(OTUB73leaf,OTUCML103leaf,OTUB73CML103leaf))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73CML103leaf)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73CML103leaf<-datarfullB73CML103leaf[,-which(readscounts<readthreshold)]
}else{
  datarB73CML103leaf<-datarfullB73CML103leaf
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73CML103leaf = otu_table(datarB73CML103leaf, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73CML103leaf[,1])){
  temp<-strsplit(rownames(datarB73CML103leaf)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73CML103leaf)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73CML103leaf<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73CML103leaf_physeq = phyloseq(OTUB73CML103leaf, TAXB73CML103leaf)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73CML103leaf_physeq)


#Get species/control-type information
B73CML103leafno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73leaf){
    B73CML103leafno<-c(B73CML103leafno,1)
  }
  else if (sample_name[k] %in% pureCML103leaf){
    B73CML103leafno<-c(B73CML103leafno,3)
  }
  else if (sample_name[k] %in% hybridB73CML103leaf){
    B73CML103leafno<-c(B73CML103leafno,2)
  }
  else{print('warning undefined maize type')}
}

B73CML103leaf_bs<-FourHbootstrap(B73CML103leaf_physeq,B73CML103leafno,0.5,500,10)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo18W Leaf 16S Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data<-read.csv(paste0('Maize/genus_bacteria.csv'), row.names=1,header= TRUE)
plantdata<-read.table('Maize/planting_design2017.txt',header=TRUE)
microbedata<-read.table('Maize/DNA_HYB.txt',header=TRUE)

pureMo18W <-plantdata$PlantID[intersect(which(plantdata$F_parent=='Mo18W'), which(plantdata$M_parent=='Mo18W'))]
hybridB73Mo18W<-plantdata$PlantID[intersect(which(plantdata$F_parent=='B73'), which(plantdata$M_parent=='Mo18W'))]

pureMo18Wmicrobe<-which(microbedata$PlantID %in% pureMo18W)
hybridB73Mo18Wmicrobe<-which(microbedata$PlantID %in% hybridB73Mo18W)

plant_part<-which(microbedata$Type=='leaf')

pureMo18Wleaf<-microbedata$SampleID[intersect(pureMo18Wmicrobe,plant_part)]
hybridB73Mo18Wleaf<-microbedata$SampleID[intersect(hybridB73Mo18Wmicrobe,plant_part)]

data<-t(data)
matrixpureB73leaf<-which(rownames(data) %in% pureB73leaf)
matrixpureMo18Wleaf<-which(rownames(data) %in% pureMo18Wleaf)
matrixhybridB73Mo18Wleaf<-which(rownames(data) %in% hybridB73Mo18Wleaf)

OTUB73leaf<-data[matrixpureB73leaf,]
OTUMo18Wleaf<-data[matrixpureMo18Wleaf,]
OTUB73Mo18Wleaf<-data[matrixhybridB73Mo18Wleaf,]

datarfullB73Mo18Wleaf<-t(rbind(OTUB73leaf,OTUMo18Wleaf,OTUB73Mo18Wleaf))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo18Wleaf)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo18Wleaf<-datarfullB73Mo18Wleaf[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo18Wleaf<-datarfullB73Mo18Wleaf
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo18Wleaf = otu_table(datarB73Mo18Wleaf, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo18Wleaf[,1])){
  temp<-strsplit(rownames(datarB73Mo18Wleaf)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo18Wleaf)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo18Wleaf<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo18Wleaf_physeq = phyloseq(OTUB73Mo18Wleaf, TAXB73Mo18Wleaf)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo18Wleaf_physeq)


#Get species/control-type information
B73Mo18Wleafno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73leaf){
    B73Mo18Wleafno<-c(B73Mo18Wleafno,1)
  }
  else if (sample_name[k] %in% pureMo18Wleaf){
    B73Mo18Wleafno<-c(B73Mo18Wleafno,3)
  }
  else if (sample_name[k] %in% hybridB73Mo18Wleaf){
    B73Mo18Wleafno<-c(B73Mo18Wleafno,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo18Wleaf_bs<-FourHbootstrap(B73Mo18Wleaf_physeq,B73Mo18Wleafno,0.5,500,10)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xCML103 Rhizosphere ITS Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data<-read.csv(paste0('Maize/fungi.csv'), row.names=1,header= TRUE)

plant_part<-which(microbedata$Type=='rhizosphere')

pureCML103rhizosphere<-microbedata$SampleID[intersect(pureCML103microbe,plant_part)]
hybridB73CML103rhizosphere<-microbedata$SampleID[intersect(hybridB73CML103microbe,plant_part)]

data<-t(data)
matrixpureB73rhizosphere<-which(rownames(data) %in% pureB73rhizosphere)
matrixpureCML103rhizosphere<-which(rownames(data) %in% pureCML103rhizosphere)
matrixhybridB73CML103rhizosphere<-which(rownames(data) %in% hybridB73CML103rhizosphere)

OTUB73rhizosphere<-data[matrixpureB73rhizosphere,]
OTUCML103rhizosphere<-data[matrixpureCML103rhizosphere,]
OTUB73CML103rhizosphere<-data[matrixhybridB73CML103rhizosphere,]

datarfullB73CML103rhizosphere<-t(rbind(OTUB73rhizosphere,OTUCML103rhizosphere,OTUB73CML103rhizosphere))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73CML103rhizosphere)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73CML103rhizosphere<-datarfullB73CML103rhizosphere[,-which(readscounts<readthreshold)]
}else{
  datarB73CML103rhizosphere<-datarfullB73CML103rhizosphere
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73CML103rhizosphere = otu_table(datarB73CML103rhizosphere, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73CML103rhizosphere[,1])){
  temp<-strsplit(rownames(datarB73CML103rhizosphere)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73CML103rhizosphere)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73CML103rhizosphere<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73CML103rhizosphere_physeq = phyloseq(OTUB73CML103rhizosphere, TAXB73CML103rhizosphere)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73CML103rhizosphere_physeq)


#Get species/control-type information
B73CML103rhizosphereno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73rhizosphere){
    B73CML103rhizosphereno<-c(B73CML103rhizosphereno,1)
  }
  else if (sample_name[k] %in% pureCML103rhizosphere){
    B73CML103rhizosphereno<-c(B73CML103rhizosphereno,3)
  }
  else if (sample_name[k] %in% hybridB73CML103rhizosphere){
    B73CML103rhizosphereno<-c(B73CML103rhizosphereno,2)
  }
  else{print('warning undefined maize type')}
}

B73CML103rhizosphere_bs<-FourHbootstrap(B73CML103rhizosphere_physeq,B73CML103rhizosphereno,0.5,500,10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Maize B73xMo18W Rhizosphere ITS Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pureMo18Wrhizosphere<-microbedata$SampleID[intersect(pureMo18Wmicrobe,plant_part)]
hybridB73Mo18Wrhizosphere<-microbedata$SampleID[intersect(hybridB73Mo18Wmicrobe,plant_part)]

matrixpureB73rhizosphere<-which(rownames(data) %in% pureB73rhizosphere)
matrixpureMo18Wrhizosphere<-which(rownames(data) %in% pureMo18Wrhizosphere)
matrixhybridB73Mo18Wrhizosphere<-which(rownames(data) %in% hybridB73Mo18Wrhizosphere)

OTUB73rhizosphere<-data[matrixpureB73rhizosphere,]
OTUMo18Wrhizosphere<-data[matrixpureMo18Wrhizosphere,]
OTUB73Mo18Wrhizosphere<-data[matrixhybridB73Mo18Wrhizosphere,]

datarfullB73Mo18Wrhizosphere<-t(rbind(OTUB73rhizosphere,OTUMo18Wrhizosphere,OTUB73Mo18Wrhizosphere))


#Calculate the reads in each sample
readscounts<-colSums(datarfullB73Mo18Wrhizosphere)

#Remove all samples below the threshold number of reads (sample failed to sequence, this number should be much smaller when environmental reads are removed)
if (length(which(readscounts<readthreshold))){
  datarB73Mo18Wrhizosphere<-datarfullB73Mo18Wrhizosphere[,-which(readscounts<readthreshold)]
}else{
  datarB73Mo18Wrhizosphere<-datarfullB73Mo18Wrhizosphere
}

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUB73Mo18Wrhizosphere = otu_table(datarB73Mo18Wrhizosphere, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datarB73Mo18Wrhizosphere[,1])){
  temp<-strsplit(rownames(datarB73Mo18Wrhizosphere)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUB73Mo18Wrhizosphere)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXB73Mo18Wrhizosphere<-tax_table(taxmat)



#Create a phyloseq object by combining the OTU table with the taxonomy object
B73Mo18Wrhizosphere_physeq = phyloseq(OTUB73Mo18Wrhizosphere, TAXB73Mo18Wrhizosphere)

#Sample identity (e.g., lizard/site)
sample_name<-sample_names(B73Mo18Wrhizosphere_physeq)


#Get species/control-type information
B73Mo18Wrhizosphereno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% pureB73rhizosphere){
    B73Mo18Wrhizosphereno<-c(B73Mo18Wrhizosphereno,1)
  }
  else if (sample_name[k] %in% pureMo18Wrhizosphere){
    B73Mo18Wrhizosphereno<-c(B73Mo18Wrhizosphereno,3)
  }
  else if (sample_name[k] %in% hybridB73Mo18Wrhizosphere){
    B73Mo18Wrhizosphereno<-c(B73Mo18Wrhizosphereno,2)
  }
  else{print('warning undefined maize type')}
}

B73Mo18Wrhizosphere_bs<-FourHbootstrap(B73Mo18Wrhizosphere_physeq,B73Mo18Wrhizosphereno,0.5,500,10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 2d

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FourHquaternary(B73Mo17leaf_bs,col='lawngreen')
FourHquaternary(B73Mo17rhizosphere_bs,col='burlywood3',add=TRUE)
FourHquaternary(B73Mo17leaffungi_bs,col='darkgreen',add=TRUE)
FourHquaternary(B73Mo17rhizospherefungi_bs,col='brown4',add=TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 2e

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHquaternary(B73Mo17leaf_bs,col='lawngreen')
FourHquaternary(B73Mo17rhizosphere_bs,col='burlywood3',add=TRUE)
FourHnullplane(B73Mo17leaf_bs,col='lawngreen')
FourHnullplane(B73Mo17rhizosphere_bs,col='burlywood3')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 2f

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FourHquaternary(B73Mo17leaf_bs,col='red',size_bootstrap=0.5)
FourHquaternary(B73CML103leaf_bs,col='gold',size_bootstrap=0.5, add=TRUE)
FourHquaternary(B73Mo18Wleaf_bs,col='blue',size_bootstrap=0.5, add=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Columns 1,2,4,5 of Table 1 and Column 1 of Table 2 for Maize lines

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHcentroid(B73Mo17leaf_bs)
FourHcentroid(B73Mo17rhizosphere_bs)
FourHcentroid(B73Mo17leaffungi_bs)
FourHcentroid(B73Mo17rhizospherefungi_bs)
FourHcentroid(B73CML103leaf_bs)
FourHcentroid(B73Mo18Wleaf_bs)
FourHcentroid(B73CML103rhizosphere_bs)
FourHcentroid(B73Mo18Wrhizosphere_bs)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Columns 2,3 of Table 2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHnullplaneD(B73Mo17leaf_bs)
FourHnullplaneD(B73Mo17rhizosphere_bs)
FourHnullplaneD(B73Mo17leaffungi_bs)
FourHnullplaneD(B73Mo17rhizospherefungi_bs)
FourHnullplaneD(B73CML103leaf_bs)
FourHnullplaneD(B73Mo18Wleaf_bs)
FourHnullplaneD(B73CML103rhizosphere_bs)
FourHnullplaneD(B73Mo18Wrhizosphere_bs)

