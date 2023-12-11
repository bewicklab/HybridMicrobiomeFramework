library("phyloseq")
library('HybridMicrobiomes')
library('stringr')
library("microViz")

#Note that numbers will differ slightly from Camper et al., due to variation in random sampling, both during matrix
#rarefaction as well as bootstrapping

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Cicada Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in OTU table
cicadamicrobiome <- read.csv("Cicadas/table_level6_cicadas.csv",row.names=1, header=TRUE)
#File relating sample ID in OTU table to SRR sequencing run number
cicadaconversion <- read.csv("Cicadas/manifestfile.csv")
#File relating SRR sequencing run number to JK number
metadata<-read.csv('Cicadas/SRRmapping.csv')
#File relating JK number to hybrid status (parent 1 ('M'), parent 2 ('T'), hybrid ('H'))
cicadaconversion2 <- read.csv("Cicadas/cicadatype.csv")

#SRR sequencing run number from manifestfile.csv
RuntoSample<-cicadaconversion$absolute.filepath
#sample ID from the manifestfile.csv
SampletoRun<-cicadaconversion$sampleid

#SRR sequencing run number from SRRmapping.csv
RuntoJK<-metadata$Run
#JK number from SRRmapping.csv
JKtoRun<-metadata$LibraryName
#JK number from cicadatype.csv
JKtoType<-cicadaconversion2$Sample
#hybrid status from cicadatype.csv
TypetoJK<-cicadaconversion2$Type
#sample ID from the OTU table
FocusSample<-colnames(cicadamicrobiome)

#for each sample ID...
typesample<-c()
cicadano<-c()
for (k in 1:length(FocusSample)){
  #find the sample ID in the list from manifestfile.csv
  hit<-which(SampletoRun==FocusSample[k])
  #find the corresponding SRR sequencing run number from manifestfile.csv
  runfocus<-RuntoSample[hit[1]]
  #find the SRR sequencing run number in the list from SRRmapping.csv
  hit2<-which(RuntoJK==runfocus)
  #find the corresponding JK number from SRRmapping.csv
  jkfocus<-JKtoRun[hit2]
  #find the JK number in the list from cicadatype.csv
  hit3<-which(JKtoType==jkfocus)
  #find the corresponding hybrid type from cicadatype.csv
  typesample<-c(typesample,TypetoJK[hit3])
  #create a numeric vector where the 'M' cicadas are assigned 1, the 'T' cicadas are assigned 3, and the hybrid 'H' cicadas are assigned 2
  if (TypetoJK[hit3]=='M'){
    cicadano<-c(cicadano,1)
  }
  else if (TypetoJK[hit3]=='H'){
    cicadano<-c(cicadano,2)
  }
  else if (TypetoJK[hit3]=='T'){
    cicadano<-c(cicadano,3)
  }
  else{
    print('warning undefined cicada type')
  }
}


#Convert dataframe into an otu table (for creating a phyloseq object)
OTUcicada = otu_table(cicadamicrobiome, taxa_are_rows = TRUE)

#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(cicadamicrobiome[,1])){
  temp<-strsplit(rownames(cicadamicrobiome)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(cicadamicrobiome)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXcicada<-tax_table(taxmat)

#Create a phyloseq object by combining the OTU table with the taxonomy object
cicada_physeq = phyloseq(OTUcicada, TAXcicada)

cicada_bs<-FourHbootstrapA(cicada_physeq,cicadano,0.5,500,7)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Lizard Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in OTU table
lizardmicrobiome <- read.csv("Lizards/Lizards.csv",row.names=1, header=TRUE)

#Find the sample names in the csv file
csv_label<-colnames(lizardmicrobiome)

lizardno<-c()
for (k in 1:length(csv_label)){
  if (substr(csv_label[k],1,1)=='I'){
    lizardno<-c(lizardno,1)
  }
  else if (substr(csv_label[k],1,1)=='N'){
    lizardno<-c(lizardno,2)
  }
  else if (substr(csv_label[k],1,1)=='M'){
    lizardno<-c(lizardno,3)
  }
  else{
    print('warning undefined lizard type')
  }
}

#Convert dataframe into an otu table (for creating a phyloseq object)
OTUlizard = otu_table(lizardmicrobiome, taxa_are_rows = TRUE)

#Create a taxonomy table (the zymo biom file doesn't seem to include this...)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
#Extract the taxonomic name at each taxonomic level
for (k in 1:length(lizardmicrobiome[,1])){
  temp<-strsplit(rownames(lizardmicrobiome)[k],split=';')
  domain<-c(domain,substr(temp[[1]][1],4,50))
  if (temp[[1]][2]=="p__NA"){temp[[1]][2]<-paste0('unclassified_',substr(temp[[1]][1],4,50))}
  else {temp[[1]][2]<-substr(temp[[1]][2],4,50)}
  phylum<-c(phylum,temp[[1]][2])
  if (temp[[1]][3]=="c__NA"){temp[[1]][3]<-paste0('unclassified_',str_replace(temp[[1]][2],'unclassified_',''))}
  else {temp[[1]][3]<-substr(temp[[1]][3],4,50)}
  class<-c(class,temp[[1]][3])
  if (temp[[1]][4]=="o__NA"){temp[[1]][4]<-paste0('unclassified_',str_replace(temp[[1]][3],'unclassified_',''))}
  else {temp[[1]][4]<-substr(temp[[1]][4],4,50)}
  order<-c(order,temp[[1]][4])
  if (temp[[1]][5]=="f__NA"){temp[[1]][5]<-paste0('unclassified_',str_replace(temp[[1]][4],'unclassified_',''))}
  else {temp[[1]][5]<-substr(temp[[1]][5],4,50)}
  family<-c(family,temp[[1]][5])
  if (temp[[1]][6]=="g__NA"){temp[[1]][6]<-paste0('unclassified_',str_replace(temp[[1]][5],'unclassified_',''))}
  else {temp[[1]][6]<-substr(temp[[1]][6],4,50)}
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}
#bind the taxonomic level names together into a table/matrix
taxmat = cbind(domain,phylum,class,order,family,genus,species)
#The rownames are the rownames from the lizardmicrobiome table
rownames(taxmat) <- rownames(lizardmicrobiome)
#The column names are the taxonomic levels
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
#Convert the taxonomy table to the format required for a phyloseq object
TAXlizard<-tax_table(taxmat)

#Create a phyloseq object by combining the OTU table with the taxonomy object
lizard_physeq = phyloseq(OTUlizard, TAXlizard)

lizard_bs<-FourHbootstrapA(lizard_physeq,lizardno,0.5,500,7)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Make Woodrat Phyloseq Object and Numeric Hybrid-Type Vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

woodratmicrobiome <- read.csv("Woodrats/table_level6.csv",row.names=1, header=TRUE)
woodratconversion <- read.csv("Woodrats/manifest_conversion.csv")
woodratconversion2 <- read.csv("Woodrats/runtable.csv")
metadata<-read.csv('Woodrats/woodratmetadata.csv')

#SRR sequencing run number from runtable.csv
Run<-woodratconversion2$Run
#S0 sequencing number from runtable.csv
Library<-woodratconversion2$Library.Name
#sample ID from manifest_conversion.csv
mylist1<-woodratconversion[,1]
#SRR sequencing run number from manifest_conversion.csv
mylist2<-woodratconversion[,2]
#S0 sequencing number from woodratmetadata.csv
Sample<-metadata$Sample
#hybrid status from woodratmetadata.csv
genotype<-metadata$Genotype

#For every sample in the order it appears in woodratmetadata.csv...
mysamplename<-c()
for (k in 1:length(Sample)){
  #Find the S0 number from runtable.csv
  hit<-which(Library==Sample[k])
  #Use the SRR number from runtable.csv to find the corresponding entry in manifest_conversion.csv
  hit2<-which(mylist2==Run[hit])
  #Find the corresponding sample name in manifest_conversion.csv
  mysamplename<-c(mysamplename,mylist1[hit2[1]])
}

#Make lists of the sample names corresponding to each hybrid/parent classification
F1s<-mysamplename[which(genotype=='F1')]
Nlbcs<-mysamplename[which(genotype=='N. lepida-bc')]
Nbbcs<-mysamplename[which(genotype=='N. bryanti-bc')]
Nls<-mysamplename[which(genotype=='N. lepida')]
Nbs<-mysamplename[which(genotype=='N. bryanti')]

#Make lists of the columns in the OTU table corresponding to each hybrid/parent classification
F1is<-which(colnames(woodratmicrobiome) %in% F1s)
Nlbcis<-which(colnames(woodratmicrobiome) %in% Nlbcs)
Nbbcis<-which(colnames(woodratmicrobiome) %in% Nbbcs)
Nlis<-which(colnames(woodratmicrobiome) %in% Nls)
Nbis<-which(colnames(woodratmicrobiome) %in% Nbs)

#Make a matrix of the different microbiome samples in order of hybrid/parent classification
datar<-cbind(woodratmicrobiome[,Nlis],woodratmicrobiome[,Nlbcis],woodratmicrobiome[,F1is],woodratmicrobiome[,Nbbcis],woodratmicrobiome[,Nbis])

#Convert dataframe into an otu table (for phyloseq - mainly used for rarifying)
OTUwoodrat = otu_table(datar, taxa_are_rows = TRUE)


#Create a taxonomy table object (for phyloseq - again mainly used for rarifying)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
for (k in 1:length(datar[,1])){
  temp<-strsplit(rownames(datar)[k],split=';')
  domain<-c(domain,temp[[1]][1])
  phylum<-c(phylum,temp[[1]][2])
  class<-c(class,temp[[1]][3])
  order<-c(order,temp[[1]][4])
  family<-c(family,temp[[1]][5])
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


taxmat = cbind(domain,phylum,class,order,family,genus,species)
rownames(taxmat) <- rownames(OTUwoodrat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
TAXwoodrat<-tax_table(taxmat)

#Create a phyloseq object by combining the OTU table with the taxonomy object
woodrat_physeq = phyloseq(OTUwoodrat, TAXwoodrat)

sample_name<-sample_names(woodrat_physeq)

#Get species
woodratno<-c()
for (k in 1:length(sample_name)){
  if (sample_name[k] %in% Nls){
    woodratno<-c(woodratno,1)
  }
  else if (sample_name[k] %in% Nlbcs){
    woodratno<-c(woodratno,6)
  }
  else if (sample_name[k] %in% F1s){
    woodratno<-c(woodratno,2)
  }
  else if (sample_name[k] %in% Nbbcs){
    woodratno<-c(woodratno,6)
  }
  else if (sample_name[k] %in% Nbs){
    woodratno<-c(woodratno,3)
  }
  else {
    woodratno<-c(woodratno,6)
  }
}

woodrat_bs<-FourHbootstrapA(woodrat_physeq,woodratno,0.5,500,7)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 3a

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHquaternary(cicada_bs,col='black',)
FourHquaternary(lizard_bs,col='green',add=TRUE)
FourHquaternary(woodrat_bs,col='brown',add=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 3b

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHquaternary(woodrat_bs,col='brown')
FourHnullplane(woodrat_bs,col='brown')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Plot Figure 3c

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHquaternary(lizard_bs,col='green')
FourHnullplane(lizard_bs,col='green')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Columns 1,2,4,5 of Table 3 and Column 1 of Table 4 for Maize lines

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHcentroid(woodrat_bs)
FourHcentroid(lizard_bs)
FourHcentroid(cicada_bs)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                            Columns 2,3 of Table 4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FourHnullplaneD(woodrat_bs)
FourHnullplaneD(lizard_bs)
FourHnullplaneD(cicada_bs)
