# A Global Analysis of multinational corporations' role in environmental conflicts #
# Marcel Llavero-Pasquina - marcelllaveropasquina@gmail.com #

#### FUNCTIONS ####

add.n <- function(factor){
  factor <- as.factor(factor)
  if(length(grep("n = ", levels(factor))) >1){}else{
    for(i in 1:length(levels(factor))){
      levels(factor)[i] <- paste0(levels(factor)[i], " (n = ", table(factor)[i], ")")
    }
  }
  return(factor)
}

add.n.cleanpoints <- function(factor){
  #Substituting point in names
  colnames(factor) <- gsub("\\.+"," ",colnames(factor))
  #Adding total n
  for(i in 2:length(colnames(factor))){
    colnames(factor)[i] <- paste0(colnames(factor)[i]," (n = ",sum(factor[,i]),")")
  }
  return(factor)
}

#Generates a table with Chi Square test between levels of a given factor for each of the categories of another factor
chi.table<- function(Xfactor, Yfactor){
  #Xfactor:
  #The independent variable to compare in the Chi.square (ie. cases$Multinationals or cases$State)
  #Yfactor:
  #The dependent variable to compare (ie. cases$Project.Status)
  Xfactor <- as.factor(Xfactor)
  if(is.list(Yfactor)){
    Yfactor$compile <- rep(NA,length.out = nrow(Yfactor)) 
    for(i in 1:nrow(Yfactor)){
      compilation <- NULL
      for(j in 2:(which(colnames(Yfactor) %in% "Multinationals")-1)){
        if(Yfactor[i,j] == 1|Yfactor[i,j] == "V"){
          compilation <- c(compilation,colnames(Yfactor)[j]) 
        }else{}
      }
      Yfactor$compile[i] <- paste(compilation, collapse = "\n")
    }
    Ytrans <- Yfactor[which(Xfactor == levels(Xfactor)[1]),]
    Ytrans <- as.data.frame(table(unlist(strsplit(Ytrans$compile,"\n"))))
    heading <- c("factor",levels(Xfactor)[1])
    for(i in 2:nlevels(Xfactor)){
      Ytrans1 <- Yfactor$compile[which(Xfactor == levels(Xfactor)[i])]
      Ytrans <- cbind.data.frame(Ytrans,as.data.frame(table(unlist(strsplit(Ytrans1,"\n"))))[,-1])
      heading <- c(heading, levels(Xfactor)[i])
    }
    #Add total n for each outcome
    levels(Ytrans$Var1) <- paste0(levels(Ytrans$Var1)," (n = ",table(unlist(strsplit(Yfactor$compile,"\n"))),")")
    colnames(Ytrans) <- heading
    chitable <- Ytrans
    
  }else{
    if(length(grep("\n",Yfactor)) > 0){
      Ytrans <- Yfactor[which(Xfactor == levels(Xfactor)[1])]
      Ytrans <- as.data.frame(table(unlist(strsplit(Ytrans,"\n"))))
      Ytrans$XFactor <- levels(Xfactor)[1] 
      heading <- c("factor",paste(levels(Xfactor)[1]))
      for(i in 2:nlevels(Xfactor)){
        Ytrans1 <- Yfactor[which(Xfactor == levels(Xfactor)[i])]
        Ytrans1 <- as.data.frame(table(unlist(strsplit(Ytrans1,"\n"))))
        Ytrans1$XFactor <- levels(Xfactor)[i]
        Ytrans <- rbind(Ytrans,Ytrans1)
        heading <- c(heading,paste(levels(Xfactor)[i]))
      }
      Xfac <- unique(Ytrans$Var1)
      data <- NULL
      XFactor <- NULL
      for(i in 1:length(Xfac)){
        num<-NULL
        for(j in 1:nlevels(Xfactor)){
          if(length(which((Ytrans$XFactor %in% levels(Xfactor)[j]) & (Ytrans$Var1 %in% levels(Xfac)[i])))>0){
            num1 <- Ytrans$Freq[which((Ytrans$XFactor %in% levels(Xfactor)[j]) & (Ytrans$Var1 %in% levels(Xfac)[i]))]
          } else {num1 <- 0}
          num <- c(num,num1)
          XFactor <- c(XFactor, levels(Xfactor)[j])
        }
        data <- rbind(data,num)
      }
      YFactor <- Xfac
      Ytrans <- cbind.data.frame(YFactor,data)
      #Add total n for each outcome
      for(i in 1:nlevels(Ytrans$YFactor)){
        levels(Ytrans$YFactor)[i] <- paste0(levels(Ytrans$YFactor)[i]," (n = ",sum(Ytrans[i,2:ncol(Ytrans)]),")")
      }
      colnames(Ytrans) <- heading
      chitable <- Ytrans
    }else{Yfactor <- as.factor(Yfactor)
    Yfactor <- add.n(Yfactor)
    chitable <- as.data.frame(levels(Yfactor))
    heading <- "factor"
    #deparse(substitute(Xfactor)) 
    for(i in 1:nlevels(Xfactor)){
      chitable <- cbind.data.frame(chitable,as.integer(table(Yfactor[which(Xfactor %in% levels(Xfactor)[i])])))
      heading <- c(heading, levels(Xfactor)[i])
    }
    colnames(chitable) <- heading}}
  
  chitable$Ratio <- as.numeric(rep(NA,length.out = nrow(chitable)))
  
  for(i in 1:nrow(chitable)){
    chitable$Ratio[i] <- chitable[i,2]/sum(chitable[i,3:ncol(chitable)],na.rm =T)
  }
  
  chisq <- chisq.test(chitable[,c(-1,-ncol(chitable))])
  chitable <- cbind(chitable,chisq$residuals)
  
  residuals <- rep(NA,length.out = nrow(chitable))
  for(i in 1:nlevels(Xfactor)){
    res<-NULL
    for(j in 1:nrow(chitable)){
      res <- c(res,as.numeric(table(Xfactor)[i]) - chitable[j,i+1])
    }
    residuals <- cbind.data.frame(residuals,res)
  }
  residuals <- residuals[,-1]
  
  chitable$pvalue <- rep(NA,length.out = nrow(chitable))
  
  for(i in 1:nrow(chitable)){
    Pretest <- NULL
    for(j in 1:nlevels(Xfactor)){
      Sample <- as.numeric(c(chitable[i,j+1],residuals[i,j]))
      Pretest <- rbind.data.frame(Pretest,Sample)
    }
    Test <- chisq.test(Pretest)
    chitable$pvalue[i] <- Test$p.value
  }
  expected <- chisq$expected
  colnames(expected) <- gsub("\\b$"," \\(expected\\)",colnames(expected))
  chitable <- cbind(chitable,expected)
  relative_abundance <- chitable[,1:nlevels(Xfactor)+1]/chitable[,(ncol(chitable)-nlevels(Xfactor)+1):ncol(chitable)]
  colnames(relative_abundance) <- gsub("\\b$"," \\(relative_abundance\\)",colnames(relative_abundance))
  chitable <- cbind(chitable,relative_abundance)
  chitable <- chitable[order(chitable$Ratio, decreasing = TRUE), ]
  return(chitable)
}

#Plots the result of the chi.table function
chi.plot <- function(chitable,Xfactors, style = "proportion"){
  Xfactors <- as.factor(Xfactors)
  total <- length(Xfactors) 
  #Xfactors <- Xfactors[-which(Xfactors %in% levels(Xfactors)[nlevels(Xfactors)])]
  chitable$Yfactor <- factor(chitable$factor, levels = c(as.character(chitable$factor)))
  data <- NULL
  Xfactor <- NULL
  yintercept <- NULL
  a <- 0
  #n <- (ncol(chitable)-4)/2
  n <- nlevels(Xfactors)
  if(style %in% c("bars","points"))
  {chitable <- chitable[order(chitable$pvalue, decreasing = TRUE), ]
  for(i in 1:n){
    data <- c(data,chitable[,(ncol(chitable)-n-1+i)])
    Xfactor <- c(Xfactor,rep(colnames(chitable)[i+1],nlevels(chitable$Yfactor)))
    a <- a+table(Xfactors)[i]
  }
  Yfactor <- rep(chitable$Yfactor,n)
  chitable_data <- cbind.data.frame(Yfactor,data)
  chitable_data$Xfactor <- Xfactor
  significant <- NULL
  for(i in 1:nrow(chitable)){
    significant <- c(significant,if(chitable$pvalue[i] > 0.05) {"no"} else {"yes"})
  }
  chitable_data$significant <- significant
  chitable_data$Xfactor <- as.factor(chitable_data$Xfactor)
  for(i in 1:nrow(chitable_data)){
    if(chitable_data$data[i] == 0){chitable_data$data[i] <-NA}
    else{
      chitable_data$data[i] <- log2(chitable_data$data[i])
    }
  }
  if(style == "bars")
  {chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, fill=Xfactor)) +
    geom_col(position=position_dodge()) +
    labs(x = NULL, y = NULL)  +
    theme(legend.position = "bottom") +
    scale_fill_manual(name = paste(Xfactors),values=c("#1e818b","#e2c537","#3BB449","#AF3456","#ED2224","#734C20","#111111","#9B989B","#BE85BA"))+
    coord_flip()+
    #gghighlight(chitable_data$significant %in% "yes",
    #unhighlighted_params =  list(fill = NULL, alpha = 0.5))+
    geom_hline(yintercept = 0, linewidth=0.5)}
  if(style == "points")
  {a <- max(na.omit(abs(chitable_data$data)))
  chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, colour=Xfactor)) +
    geom_point(size = 6) +
    labs(x = NULL, y = NULL)  +
    theme(legend.position = "bottom") +
    scale_colour_manual(values=c("#1e818b","#e2c537","#3BB449","#AF3456","#ED2224","#734C20","#111111","#9B989B","#BE85BA"))+
    coord_flip()+
    gghighlight(chitable_data$significant %in% "yes",
                unhighlighted_params =  list(colour = NULL, alpha = 0.3))+
    geom_hline(yintercept = 0, linewidth=0.5)+
    ylim(c(-a,a))
  #geom_text(hjust=0.5, vjust=0.5, colour="#ffffff", fontface = "bold")
  }
  chiplot
  }
  else
  {for(i in 1:n){
    data <- c(data,chitable[,i+1])
    Xfactor <- c(Xfactor,rep(colnames(chitable)[i+1],nlevels(chitable$Yfactor)))
    a <- a+table(Xfactors)[i]
    yintercept <- c(yintercept,a/total)
  }
    yintercept <- 1-yintercept
    Yfactor <- rep(chitable$Yfactor,n)
    chitable_data <- cbind.data.frame(Yfactor,data)
    chitable_data$Xfactor <- Xfactor
    significant <- NULL
    for(i in 1:nrow(chitable)){
      significant <- c(significant,if(chitable$pvalue[i] > 0.05) {"no"} else {"yes"})
    }
    chitable_data$significant <- significant
    
    chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, fill=Xfactor)) +
      geom_col(position="fill",) +
      labs(x = NULL, y = NULL)  +
      theme(legend.position = "bottom") +
      scale_fill_manual(name = paste(Xfactors),values=c("#1e818b","#e2c537","#F6EB13","#7A752F","#734C20","#111111","#9B989B","#BE85BA","#3BB449","#ED2224","#AF3456"))+
      coord_flip()+
      gghighlight(chitable_data$significant %in% "yes",
                  unhighlighted_params =  list(fill = NULL, alpha = 0.5))+
      geom_hline(yintercept = yintercept, linewidth=0.5)}
  
  return(chiplot)
}

save.plot <- function(chiplot){
  ggsave(chiplot, filename = paste0(deparse(substitute(chiplot)),".svg"), width = 7, height = (1+nlevels(chiplot$data$Yfactor)/5), dpi = 300, units = "in")
}

#### LIBRARIES ####
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(gghighlight)
library(data.table)
library(remotes)
library(ggsankey)
library(svglite)

#### 1. Import and prepare data ####

setwd("C:/Users/Usuario/Desktop/Marcel/EJAtlas/Writings/Political Ecology Theory of the Firm/EJAtlas Company Analysis/Data")
cases<-read.csv("cases.csv", sep = ",")
companys<-read.csv("companys.csv", sep = ",")
company_conflicts<-read.csv("company_conflicts.csv", sep = ",")
conflict_companys<-read.csv("conflict_companys.csv", sep = ",")
mobilizinggroups <- read.csv("mobilizinggroups.csv", sep = ",")
mobilizingforms <- read.csv("mobilizingforms.csv", sep = ",")
envimpacts <- read.csv("envimpacts.csv", sep = ",")
hltimpacts <- read.csv("hltimpacts.csv", sep = ",")
secimpacts <- read.csv("secimpacts.csv", sep = ",")
Regions <- read.csv("Regions.csv")

Ranking <- companys

#Correct a couple of identified spelling mistakes
Ranking$name[which(Ranking$name %in% "BPH Billiton")] <- "BHP Billiton"
Ranking$name[which(Ranking$name %in% "Xtrata Copper")] <- "Xstrata Copper"

#Clean all spaces, points and commas at beginning or end of name
Ranking$clean_name <- gsub("^[ \\.,]*|[ \\.,]*$", "", Ranking$name)

#### 2. Identify most common words in company names ie. "National" or "Mining ####

Word <- as.data.frame(table(unlist(strsplit(companys$name, "\\s+"))))
table(Word$Freq)

#### 3. Set specificity for the company matching ####
#The specificity value is the exclusion threshold for the number of 
#repetitions of a word in the companies names.
#ie. any word repeated more than 9 times in all company names will not
#be considered for matching companies names.
#The lower the value the more false positive matches, 
#the higher the value the more missed high-impact companies like Chevron or Shell.
#However, the threshold can be increased, while manually parsing high-impact companies.

Specificity <- 2
Words <- as.character(Word$Var1[which(Word$Freq > Specificity)])

#Also include all country names for exclusion
library(countrycode)
a<-countrycode::codelist[,6]
Words <- c(Words,a)
rm(a)

#Append other generic words or words that create problems later on manually
Words<- c(Words,
          "energy",
          "mining",
          "Inter",
          "One",
          "Petro",
          "Invest",
          "Hydro",
          "Trans",
          "Consult",
          "minera",
          "Anglo",
          "Corporacion",
          "Hydroelectric",
          "Electricals",
          "PLC",
          "Petroleos",
          "Petróleo",
          "Corp.",
          "Corp",
          "Inc.",
          "inc.",
          "Ltd",
          "Andhra",
          "Hindustan",
          "Tokyo",
          "Québec",
          "Saudi",
          "Nile",
          "Energias",
          "Offshore",
          "&",
          "Gran",
          "Motors",
          "Silver",
          "Harbour",
          "Communication",
          "Atómica",
          "Organization",
          "Americas",
          "Jatropha",
          "Combustibles",
          "Portfolio",
          "Diamond",
          "INC",
          "consortium",
          "INGENIEROS",
          "Pty",
          "S. A",
          "Ltd",
          "SEA",
          "Rica",
          "East",
          "Cia",
          "Geo",
          "Pte",
          "DI",
          "CPI",
          "SRI",
          "Congolaise",
          "Navy",
          "Ingeniería",
          "Sai",
          "Owned",
          "ICE",
          "Dutch",
          "Hellas",
          "Leone",
          "Mocambique",
          "Brasileiro",
          "Reserve",
          "Comission",
          "Portugal",
          "Eko",
          "Hidrelétrica",
          "Laos",
          "Fondo",
          "Greenland",
          "Española",
          "Detroit",
          "Harvest",
          "Base",
          "Singapore",
          "Bulk",
          "Fomento",
          "Contractors",
          "Imperial",
          "Austral",
          "Metalúrgica",
          "Fujian",
          "Ind",
          "Assam",
          "Noroeste",
          "Mutare",
          "Verde",
          "Bogotá",
          "Fé",
          "Titanium",
          "Sarawak",
          "Flow",
          "Corridor",
          "Sucursal",
          "Kingho",
          "Santa Rita",
          "Ship",
          "Breakers",
          "Jilin",
          "Switzerland",
          "Powergrid",
          "Arabian",
          "Facilities",
          "Venezolana",
          "Senegal",
          "Aurora",
          "Société",
          "Manuelita",
          "Asturias",
          "Man",
          "Emirates",
          "Wastewater",
          "Northwest",
          "Vega",
          "Scotia",
          "Odisha",
          "Brasileiras",
          "Peaks",
          "Alberta",
          "Defence",
          "Croatia",
          "Energetica",
          "Gaia",
          "Liberia",
          "Manila",
          "Building",
          "Pampa",
          "Lake",
          "Santander",
          "Volta",
          "Arizona",
          "Arco",
          "Luliang",
          "Barracuda",
          "Southwest",
          "Third"
          )

#### 4. Create a clean name for companies omitting repetitive words ####

#Clean parentheses as they give regex errors
Words <- gsub("\\(|\\)", "", Words)
for(i in 1:nrow(Ranking)){
  Ranking$clean_name[i] <- gsub("\\(|\\)", "", Ranking$clean_name[i])
}
rm(i)

#Remove repetitive words from company names
Ranking$clean_name <- gsub(ignore.case = TRUE, paste0("\\b",unlist(strsplit(Words[1:500], " ")),"\\b", collapse = "|"), "", Ranking$clean_name)
Ranking$clean_name <- gsub(ignore.case = TRUE, paste0("\\b",unlist(strsplit(Words[501:1000], " ")),"\\b", collapse = "|"), "", Ranking$clean_name)
Ranking$clean_name <- gsub(ignore.case = TRUE, paste0("\\b",unlist(strsplit(Words[1001:length(Words)], " ")),"\\b", collapse = "|"), "", Ranking$clean_name)

#Clean all spaces, points and commas at beginning or end of name
Ranking$clean_name <- gsub("^[ \\.,]*|[ \\.,]*$", "", Ranking$clean_name)
#Clean all one or two letter names
Ranking$clean_name <- gsub("^.$|^..$|^.\\..$", "", Ranking$clean_name)

#If all words in the name are repetitive words, keep the original name
for(i in 1:nrow(Ranking)){
  if(grepl("\\S",Ranking$clean_name[i])){} else {
    Ranking$clean_name[i] <- Ranking$name[i]
  }  
}
rm(i)

#Fix company entry without name
Ranking$clean_name[which(Ranking$clean_name %in% "")] <- "No_name"

#Clean companies with unclosed parenthesis that give errors
Ranking$clean_name[which(Ranking$id %in% 5695)] <- "Tianji"
Ranking$clean_name[which(Ranking$id %in% 3105)] <- "EDRACO . Member  Edrasis - C. Psallidas"
Ranking$clean_name[which(Ranking$id %in% 7014)] <- "China Coal Energy Group "

#Clean a few acronyms that give false positive hits
acronymDELETE <- c("EM","CPC","GR","ADM","IDC","NCC","TEAM","CMDC","CCC","ORE","AMA",
                   "UAM","Max","EMAS","(CES)","DMP","GEM","EDM","IMC","UCI","EN+","CFA",
                   "TMC","STC","OCP","CGN","AC","MAC","SMEC","BFI","WTI","USA","ENEX","BK",
                   "ADP","GCO","GCM","PPG","SPIC","HECO","SCA","CIL","SMC","NTPC","AES")
Ranking$acronym[which(Ranking$acronym %in% acronymDELETE)] <- ""
rm(acronymDELETE)

#### 5. Add manual search word ####

#Create template spreadsheet to fill manually
Manual_search <- as.data.frame(Ranking$clean_name,Ranking$id)
Manual_search$Manual_search <- as.character(rep(NA,nrow(Ranking)))
write.csv(Manual_search,"Manual_search_to_fill.csv")

#Add words in excel and reimport
#If more than one search query is needed it is separated by /
Manual_search<-read.csv("Manual_search.csv", sep = ",")

#Pass manual_search strings into Ranking data frame based on company ID
Ranking$manual_search <- as.character(rep(NA,nrow(Ranking)))
for(i in 1:nrow(Ranking)){
  if (length(Manual_search$manual_search[which(Manual_search$id %in% Ranking$id[i])]) > 0) {
    Ranking$manual_search[i] <- Manual_search$manual_search[which(Manual_search$id %in% Ranking$id[i])]
  } else {}
}
rm(i)

#Clean NA introduced which can give problems later on
Ranking$manual_search[which(is.na(Ranking$manual_search))] <- ""

#### 6. Match company names ####

Matches <- NULL
for(i in 1:nrow(Ranking)){#!!!WARNING: EXPECT A LONG RUN TIME IN A PC!!!
  Matches <- c(Matches,paste(collapse = " ",
  unique(c(fromLast = FALSE,
                     
  #Search with exact word clean_name in name
  Ranking$id[grep(paste( 
  "\\b",Ranking$clean_name[i],"\\b", 
  sep = ""), Ranking$name, ignore.case = TRUE)]
                     
  #Search with exact word clean_name without spaces or points in name
  ,Ranking$id[grep(paste(
  "\\b",gsub(" |\\.", "",Ranking$clean_name[i]),"\\b", 
  sep = ""), Ranking$name, ignore.case = TRUE)]
                     
  #Search with exact word clean_name in acronym
  ,Ranking$id[grep(paste(
  "\\b",Ranking$clean_name[i],"\\b", 
  sep = ""), Ranking$acronym, ignore.case = TRUE)]
                     
  #Search with exact word clean_name without spaces or points in acronym
  ,Ranking$id[grep(paste( 
  "\\b",gsub(" |\\.", "",Ranking$clean_name[i]),"\\b", 
  sep = ""), Ranking$acronym, ignore.case = TRUE)]
  
  #Search with exact word acronym in name
  ,if(Ranking$acronym[i] %in% c("",Words)){} else {Ranking$id[grep(paste( 
  "\\b",gsub(" |\\.", "",Ranking$acronym[i]),"\\b", 
  sep = ""), Ranking$name, ignore.case = TRUE)]}
                     
  #Search with exact word acronym in acronym
  ,if(Ranking$acronym[i] %in% c("",Words)){} else {Ranking$id[grep(paste( 
  "\\b",gsub(" |\\.", "",Ranking$acronym[i]),"\\b", 
  sep = ""), Ranking$acronym, ignore.case = TRUE)]}
                   
  #Search with exact manual_Search words (divided by "/") in name
  ,if(Ranking$manual_search[i] %in% ""){} else {
  Ranking$id[grep(paste0("\\b",unlist(strsplit(Ranking$manual_search[i], "/")),"\\b", collapse = "|"),
  Ranking$name, ignore.case = TRUE)]}
                   
  #Search with exact manual_Search words (divided by "/") in acronym
  ,if(Ranking$manual_search[i] %in% ""){} else {
  Ranking$id[grep(gsub(" ","",paste("\\b",unlist(strsplit(Ranking$manual_search[i], "/")),"\\b", collapse = "|")),
  Ranking$acronym, ignore.case = TRUE)]}
    ))))
}
rm(i)

Ranking$Matches <- Matches
#Clean matches string
Ranking$Matches <- gsub("\\b0 \\b|\\b0\\b","",Ranking$Matches)

#How many hits for each row?
Ranking$Hits<-lengths(regmatches(Ranking$Matches, gregexpr("\\b^| ", Ranking$Matches)))

rm(Word, Words, Manual_search)

#### 7. Cleaning of hits for companies with common names/acronyms ####

#Remove matches for companies with too many hits, which are false positives
#Companies with common names
Ranking$Matches[which(Ranking$Hits >13)] <- Ranking$id[which(Ranking$Hits >13)] 
#Company named "Subsidiary"
Ranking$Matches[which(Ranking$id %in% 230)] <- Ranking$id[which(Ranking$id %in% 230)] 

#Recalculate hits
Ranking$Hits<-as.integer(lengths(regmatches(Ranking$Matches, gregexpr("\\b^| ", Ranking$Matches))))
table(Ranking$Hits) #Distribution of hits
sum(Ranking$Hits) #Total number of hits

#### 8. Append matching name matrix to Ranking ####

#Add matches to companyIDs list and count the number of companies for each entry
Ranking$Num_companies <- as.numeric(rep(NA,length.out = nrow(Ranking)))
Ranking$companyIDs <- as.numeric(rep(NA,length.out = nrow(Ranking)))
for(i in 1:nrow(Ranking)){
  Ranking$companyIDs[i] <- paste(unique(unlist(strsplit(c(Ranking$id[i],Ranking$Matches[i]), " "))), collapse = " ")
  Ranking$companyIDs[i] <- gsub("  "," ",Ranking$companyIDs[i])
  Ranking$Num_companies[i] <- length(c(unlist(strsplit(Ranking$companyIDs[i], " "))))
}

#It writes the names of the companies corresponding to the matched ids. To facilitate manual quality inspection
MatchingNames <- matrix(nrow = nrow(Ranking), ncol = max(Ranking$Num_companies))
Ranking <- cbind.data.frame(Ranking,MatchingNames)
rm(MatchingNames)
i<-1
for(i in 1:nrow(Ranking)){
  Ranking[i,(ncol(Ranking)-max(Ranking$Num_companies)+1):ncol(Ranking)] <- 
    c(companys$name[which(companys$id %in% c(unlist(strsplit(Ranking$companyIDs[i], " "))))],
      rep(NA,max(Ranking$Num_companies-length(companys$name[which(companys$id %in% c(unlist(strsplit(Ranking$companyIDs[i], " "))))]))))
}

#### 9. Manual quality check ####

#Create csv to manually inspect
write.csv(Ranking, "RankingDELETE_to_fill.csv")

#Manually check 1300+ lines of corresponding names
#Substitute with DELETE any false positive hits, including joint ventures, or hits associated with common words
#Matches should only account for 
##1. Same companies with different spelling
##2. Wholy-owned subsidiaries
##3. Mergers and acquisitions
##4. Name changes
#In some cases (Rosatom, PT PLN and Newmont) look through the database for subsidiaries without the same name

#Import the manually-checked data.frame
RankingDELETE<-read.csv("RankingDELETE.csv", sep = ",")
#Select only the rows with manually DELETED hits
RankingDELETE <- RankingDELETE[grepl("DELETE",do.call(paste,RankingDELETE)),]
#390/6729 entries that contain a DELETE
#Count number of manually deleted hits
DELETE <- 0
for(i in match("X1.1",colnames(RankingDELETE)):(match("X1.1",colnames(RankingDELETE))+max(RankingDELETE$Hits))){
  DELETE <- sum(DELETE,sum(RankingDELETE[,i] %in% "DELETE"))
}

#For each row, delete DELETED company IDs from companyIDs factor
for(i in 1:nrow(Ranking)){
  #Only for those rows that have a manual DELETE
  if(length(grep(Ranking$id[i],RankingDELETE$id)) > 0){
    #Collect the names of all the companies that were hit
    nameDELETE <- paste(RankingDELETE[which(RankingDELETE$id %in% Ranking$id[i]),match("X1.1",colnames(RankingDELETE)):match("X10.1",colnames(RankingDELETE))], collapse= "/")
    #Clean up string
    nameDELETE <- gsub("/NA", "", nameDELETE)
    #Identify company IDs from manually DELETED companies that need to be removed from companyIDs list
    DELETEmatches <- c(unlist(strsplit(RankingDELETE$Matches[which(RankingDELETE$id %in% Ranking$id[i])], split = " ")))
    idDELETE <- DELETEmatches[-which(DELETEmatches %in% companys$id[which(companys$name %in% c(unlist(strsplit(nameDELETE, split = "/"))))])]
    #Remove selected IDs from companyIDs list
    Ranking$companyIDs[i] <- gsub(paste0("\\b",idDELETE,"\\b", collapse = "|"),"",Ranking$companyIDs[i])
    Ranking$companyIDs[i] <- gsub("  "," ",Ranking$companyIDs[i])
  } else {}
}
rm(i,nameDELETE,idDELETE,DELETEmatches,DELETE,RankingDELETE)

#### 10. Assign conflict IDs to each company based on the companyIDs list ####

Ranking$conflictIDs <- as.character(rep(NA, length.out = nrow(Ranking)))
Ranking$Num_cases <- as.numeric(rep(NA, length.out = nrow(Ranking)))
for(i in 1:nrow(Ranking)){
  #List all cases IDs for each company, only unique case IDs
  conflicts <- paste(unique(conflict_companys[which(conflict_companys$company_id %in% c(unlist(strsplit(Ranking$companyIDs[i], " ")))),
                                              match("conflict_id_1",colnames(conflict_companys)):ncol(conflict_companys)])
                     , collapse = " ")
  #Remove NAs,c and parenthesis introduced for duplicated cases
  conflicts <- gsub("c\\(NA\\,\\)|\\,|\\)|c\\(","",conflicts)
  conflicts <- gsub(" NA","",conflicts)
  Ranking$conflictIDs[i] <- paste(unique(unlist(strsplit(conflicts, " "))), collapse = " ")
  
  #Count the number of cases for each company
  Ranking$Num_cases[i] <- length(c(unlist(strsplit(Ranking$conflictIDs[i], " "))))
}
rm(i,conflicts)
Ranking$conflictIDs <- gsub("\\s+"," ",Ranking$conflictIDs)

#### 11. Remove company entries of duplicates, subsidiaries, mergers, acquisitions, and name changes ####

#Saving the file before deleting company entries
write.csv(Ranking,"Ranking.csv")

#Assigning DELETE to those companies already hit by another company with more conflicts #
Ranking <- Ranking[order(Ranking$Num_cases,decreasing = TRUE),]
Ranking$DELETE <- rep("",length.out = nrow(Ranking))
for(i in 1:nrow(Ranking)){
  if(length(grep(paste0("\\b",Ranking$id[i],"\\b"),Ranking$companyIDs[1:i-1])) > 0){
      Ranking$DELETE[i] <- "DELETE"
  }
}

#Salvaging those companies which are hit by companies which are going to be deleted#

RowsDELETE <- which(Ranking$DELETE %in% "DELETE")
RowsNOTDELETE <- which(!Ranking$DELETE %in% "DELETE")
for(i in 1:nrow(Ranking)){
#Select only rows that are going to be deleted
    if(Ranking$DELETE[i] %in% "DELETE"){
      #Select those rows that are hit by companies which are going to be deleted#
    if(length(grep(paste0("\\b",Ranking$id[i],"\\b"),Ranking$companyIDs[RowsDELETE[-which(RowsDELETE > i-1)]])) > 0){
      #Select against those rows that are hit by companies which are not going to be deleted#
      if(length(grep(paste0("\\b",Ranking$id[i],"\\b"),Ranking$companyIDs[RowsNOTDELETE])) > 0){} else {
        Ranking$DELETE[i] <- "" 
      }
    }
  }
}
rm(i,RowsDELETE,RowsNOTDELETE)
#Delete companies which already appear under another company
Ranking <- Ranking[-which(Ranking$DELETE %in% "DELETE"),]

#### 12. Check that we are not missing any original company id ####

companys$check <- rep("NA",length.out = nrow(companys))
#Go through all original companys IDs and checks that they are included in at least one company in the output data frame
for(i in 1:nrow(companys)){
  companys$check[i] <- if(length(grep(paste0("\\b",companys$id[i],"\\b"),Ranking$companyIDs)) > 0) {""}
  else
  {"error"}
}
#List of all companies lost in the process. If it all works, it has to be 0
recover <- companys[which(companys$check %in% "error"),]
#Recover one company missing (YPF, id: 3628)
Ranking$companyIDs[which(Ranking$id %in% 3628)] <- "3628 2863 6200 6201"
Ranking <- Ranking[-which(Ranking$id %in% c("2863","6200","6201")),]
#Rerun check to make sure there are 0 errors

#Remove some ":" introduced in conflictIDS
Ranking$conflictIDs <- gsub("\\:"," ",Ranking$conflictIDs)

rm(recover)

#### 13. Select relevant factors for further analysis and save as Summary ####

Summary <- Ranking[,match(c("id","name","slug","description","url","acronym","country","companyIDs","Num_companies","conflictIDs","Num_cases"),colnames(Ranking))]

write.csv(Summary,"Summary.csv")
rm(Ranking,Matches,Specificity,i)

#### 14. Assign category to each company ####

#The most frequent conflict category is assigned as the company category
Summary$Category <- as.numeric(rep(NA,nrow(Summary)))
for(i in 1:nrow(Summary)){
  Summary$Category[i] <-names(sort(
    table(cases$First.level.category[which(cases$Conflict.Id %in% c(unlist(strsplit(Summary$conflictIDs[i], " "))))])
    ,decreasing = T))[1]
}

#### 15. Assign level reported to each company ####

Summary$Reported <- as.numeric(rep(NA,nrow(Summary)))
for(i in 1:nrow(Summary)){
  Summary$Reported[i] <- if(Summary$Num_cases[i] >30){
    "More than 30 cases"
  } else if(Summary$Num_cases[i] >7){
    "8-30 cases"
  } else if(Summary$Num_cases[i] >4){
    "5-7 cases"
  } else if(Summary$Num_cases[i] > 1){
    "2-4 cases"
  } else if(Summary$Num_cases[i] > 0){
    "One case"
  } 
}

Summary$Reported <- as.factor(Summary$Reported)
levels(Summary$Reported)
Summary$Reported <- factor(Summary$Reported, levels = c("One case", "2-4 cases", "5-7 cases", "8-30 cases","More than 30 cases"))
sort(table(Summary$Reported), decreasing =T)

#### 16. Assign country and region of origin to each company ####

#Counts the number of conflicts for each original company entry
conflict_companys$Num_Cases <- as.numeric(rep(NA,nrow(conflict_companys)))
for(i in 1:nrow(conflict_companys)){
  conflict_companys$Num_Cases[i] <- ncol(conflict_companys[i,which(!is.na(conflict_companys[i,]))])-1
}

#Assign the country of origin most frequently reported in all company entries related cases 
#Watch out! Summary$country (in lower-case) is the original country register of the entry!!
Summary$Country <- as.numeric(rep(NA,nrow(Summary)))
for(i in 1:nrow(Summary)){
  Countries <- companys$country[which(companys$id %in% c(unlist(strsplit(Summary$companyIDs[i], " "))))]
  Num_Cases <- conflict_companys$Num_Cases[which(conflict_companys$company_id %in% c(unlist(strsplit(Summary$companyIDs[i], " "))))]
  a<-cbind.data.frame(Countries,Num_Cases)
  a<- a[which(!(a$Countries %in% "")),]
  b<-a %>% group_by(Countries) %>% summarize(sum(Num_Cases))
  Summary$Country[i] <- if(nrow(a) == 0) {"ND"} else {b$Countries[which(b$`sum(Num_Cases)` %in% max(b$`sum(Num_Cases)`))]}
}
rm(a,b,Countries,Num_Cases)
Summary$Country<-gsub("NA","ND",Summary$Country)

#Reorganise WB income levels
Regions$World..Bank <- gsub("Lower Middle", "Low", Regions$World..Bank)
Regions$World..Bank <- gsub("Upper Middle", "Middle", Regions$World..Bank)

Summary$Region <- as.numeric(rep(NA,nrow(Summary)))
for(i in 1:nrow(Summary)){
  Summary$Region[i] <- if(Summary$Country[i] %in% "") {"ND"} else {
    Regions$World..Bank[which(Regions$Country %in% Summary$Country[i])]
  }
}
#The table Regions also allows to assign other attributes like geographical region or Democracy Index
Summary$Region <- as.factor(Summary$Region)
Summary$Region <- factor(Summary$Region, levels = levels(Summary$Region)[c(2,1,4,3)])
sort(table(Summary$Region), decreasing = T)

#### 17. Identify multinational companies ####

#If a company does not have a country of origin assigned "ND", it cannot be established whether it is a multinational
Summary$Multinational <- as.numeric(rep(NA,nrow(Summary)))
for(i in 1:nrow(Summary)){
  Summary$Multinational[i] <- if(Summary$Country[i] %in% "ND") {"ND"} else {
    if("FALSE" %in% as.character(cases$Country[which(cases$Conflict.Id %in% c(unlist(strsplit(Summary$conflictIDs[i], " "))))] %in% Summary$Country[i])) 
    {"yes"} else {"no"}
  }
}

sort(table(Summary$Multinational),decreasing = T)

#Save table
write.csv(Summary,"Summarywithattributes.csv")

#### 18. Clean up Intensity and Project status codes in cases table ####

table(cases$Intensity.of.Conflict)
cases$Intensity.of.Conflict[which(cases$Intensity.of.Conflict %in% "")] <- "UNKNOWN"
cases$Intensity.of.Conflict <- as.factor(cases$Intensity.of.Conflict)
levels(cases$Intensity.of.Conflict) <- c("HIGH","LATENT","LOW","MEDIUM","UNKNOWN")
cases$Intensity.of.Conflict <- factor(cases$Intensity.of.Conflict, levels = c("UNKNOWN","LATENT","LOW","MEDIUM","HIGH"))

table(cases$Project.Status)
cases$Project.Status[which(cases$Project.Status %in% "")] <- "Unknown"
cases$Project.Status <- as.factor(cases$Project.Status)
levels(cases$Project.Status) <- c("In operation","Planned","Proposed","Stopped","Under construction","Unknown")
cases$Project.Status <- factor(cases$Project.Status, levels = c("Unknown","Proposed","Planned","Under construction","In operation","Stopped"))

#### 19. Quality filter ####

#Filter out cases without companies
length(which(cases$Company %in% ""))
#There are 563 cases without companies assigned. Delete them
mobilizingforms <- mobilizingforms[-which(cases$Company %in% ""),]
mobilizinggroups <- mobilizinggroups[-which(cases$Company %in% ""),]
envimpacts <- envimpacts[-which(cases$Company %in% ""),]
hltimpacts <- hltimpacts[-which(cases$Company %in% ""),]
secimpacts <- secimpacts[-which(cases$Company %in% ""),]
cases <- cases[-which(cases$Company %in% ""),]

#Filter out cases older than 1947
cases$startyear<- as.numeric(rep("NA",length.out = nrow(cases)))
for(i in 1:nrow(cases)){
  cases$startyear[i] <- paste(regmatches(cases$Start.Date[i], gregexpr("\\d{4}",cases$Start.Date[i])))
}
cases$startyear <- as.numeric(cases$startyear)

length(which(cases$startyear < 1947))
#There are 46 cases we delete because they are earlier than 1947. Delete them
mobilizingforms <- mobilizingforms[-which(cases$startyear < 1947),]
mobilizinggroups <- mobilizinggroups[-which(cases$startyear < 1947),]
envimpacts <- envimpacts[-which(cases$startyear < 1947),]
hltimpacts <- hltimpacts[-which(cases$startyear < 1947),]
secimpacts <- secimpacts[-which(cases$startyear < 1947),]
cases <- cases[-which(cases$startyear < 1947),]

#Filter out companies without country assigned
length(which(Summary$Country %in% "ND"))
Summary <- Summary[-which(Summary$Country %in% "ND"),]
#There are 819 companies without country defined. Most of them only have one case associated.

#Filter out cases that only have companies without country assigned involved
Delete <- NULL
for(i in 1:nrow(cases)){
  if(length(grep(paste0("\\b",cases$Conflict.Id[i],"\\b"), Summary$conflictIDs)) >0){}
  else {Delete <- c(Delete,i) }
}
length(Delete)
#There are 119 conflicts that only had a company without country assigned involved
mobilizingforms <- mobilizingforms[-Delete,]
mobilizinggroups <- mobilizinggroups[-Delete,]
envimpacts <- envimpacts[-Delete,]
hltimpacts <- hltimpacts[-Delete,]
secimpacts <- secimpacts[-Delete,]
cases <- cases[-Delete,]
rm(Delete)

#### 20. Identify cases with at least one foreign company involved ####
cases$Multinationals <- as.numeric(rep(NA,nrow(cases)))
for(i in 1:nrow(cases)){
  if(cases$Conflict.Id[i] %in% company_conflicts$conflict_id){
    Comp <- as.numeric(company_conflicts[which(company_conflicts$conflict_id %in% cases$Conflict.Id[i]),seq(from = 2, length.out = 44, by = 4)], na.rm = TRUE)
    Comp <- as.character(Comp[!is.na(Comp)])
    Ctry <- NULL
    for(j in 1:length(Comp)){
      Ctry <- c(Ctry,Summary$Country[grep(paste("\\<",Comp[j],"\\>", sep = ""),Summary$companyIDs)] == cases$Country[i])
    }
    if(FALSE %in% Ctry)
    {cases$Multinationals[i] <- "yes"} else {cases$Multinationals[i] <- "no"}
  } else {cases$Multinationals[i] <- "no"}
}
rm(i,j,Comp,Ctry)
mobilizingforms$Multinationals <- cases$Multinationals
mobilizinggroups$Multinationals <- cases$Multinationals
envimpacts$Multinationals <- cases$Multinationals
hltimpacts$Multinationals <- cases$Multinationals
secimpacts$Multinationals <- cases$Multinationals
table(cases$Multinationals)

#### 21. Assign region of conflict ####

cases$Region <- as.numeric(rep(NA,nrow(cases)))
for(i in 1:nrow(cases)){
  cases$Region[i] <- if(cases$Country[i] %in% "") {"ND"} else {
    Regions$World..Bank[which(Regions$Country %in% cases$Country[i])]
  }
}  
#The table Regions also allows to assign other attributes like geographical region or Democracy Index

cases$Region <- as.factor(cases$Region)
mobilizingforms$Region <- cases$Region
mobilizinggroups$Region <- cases$Region
envimpacts$Region <- cases$Region
hltimpacts$Region <- cases$Region
secimpacts$Region <- cases$Region

sort(table(cases$Region), decreasing = T)

cases$Region <- factor(cases$Region, levels = levels(cases$Region)[c(2,1,4,3)])

#### 22. Bar plot companies by number of cases and multinational ####
Summary$Reported <- add.n(Summary$Reported)
#Introduces a line break in reported categories to improve graphic display 
#Summary$Reported <- gsub("\\(","\n\\(",Summary$Reported)
Summary$Reported <- as.factor(Summary$Reported)

Freq_Multinational_Reported <- Summary %>% group_by(Multinational,Reported) %>% summarize(count=n())

Freq_Multinational_Reported_plot <- ggplot(Freq_Multinational_Reported, aes(x = Reported, y = count, fill = Multinational)) +
  geom_col(position="fill")+
  labs(x = NULL, y = NULL)  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#1e818b","#e2c537"))
Freq_Multinational_Reported_plot

ggsave(Freq_Multinational_Reported_plot, filename = "Freq_Multinational_Reported_plot.svg")

#### 23. Sankey diagram of countries/regions of cases vs countries/regions of company origins ####

#Make list of all company involvements in cases indicating regions of company and conflict (n= 7662)
Sankey<-NULL
Summary$Region <- as.character(Summary$Region)
for(i in 1:nrow(Summary)){
  a <- as.numeric(unlist(strsplit(Summary$conflictIDs[i]," ")))
  Countries <- as.character(cases$Region[which(cases$Conflict.Id %in% a)])
  Companies <- rep(Summary$Region[i],length(Countries))
  d <- cbind.data.frame(Countries,Companies)
  Sankey <- rbind(Sankey,d)
}
rm(a,d, Countries, Companies)

#Convert to input format for ggsankey
df <- Sankey %>%
  make_long(Countries,Companies)
df$node <- as.factor(df$node)
df$node <- factor(df$node, levels = levels(df$node)[c(3,4,1,2)])
df$next_node <- as.factor(df$next_node)
df$next_node <- factor(df$next_node, levels = levels(df$next_node)[c(3,4,1,2)])

pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      ,node.color = "white"
                      ,show.legend = FALSE
)
pl <- pl + theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl + theme(axis.title = element_blank()
                 , axis.text.y = element_blank()
                 , axis.ticks = element_blank()  
                 , panel.grid = element_blank())
pl <- pl + labs(fill = 'Nodes')

pl

#Beautify and add conflict and companies names and n's with Adobe Illustrator
ggsave("Sankey.svg",pl, device="svg")

#Summary table of company region and conflict region interaction
Counts <- Sankey %>% group_by(Countries,Companies) %>% summarize(length(Companies))
write.csv(Counts,"Sankey.csv")

rm(Sankey,df,pl,Counts)

#### 24. Chi-squared tests of conflicts with vs without foreign companies ####

### Project status ###

Status <- chi.table(cases$Multinationals,cases$Project.Status)
Status_plot <-chi.plot(Status,cases$Multinationals)

### Category ###

Category <- chi.table(cases$Multinationals,cases$First.level.category)
Category_plot <- chi.plot(Category,cases$Multinationals)

### Intensity ###

Intensity <- chi.table(cases$Multinationals,cases$Intensity.of.Conflict)
Intensity_plot <- chi.plot(Intensity,cases$Multinationals)

### Outcomes ###

Outcomes <- chi.table(cases$Multinationals,cases$Outcomes)
Outcomes_plot <- chi.plot(Outcomes,cases$Multinationals)

### Second level types ###

Types <- chi.table(cases$Multinationals,cases$Second.level.type)
Types_plot <- chi.plot(Types,cases$Multinationals)

### Commodities ###

Commodities <- chi.table(cases$Multinationals,cases$Commodity)
Commodities_plot <- chi.plot(Commodities,cases$Multinationals)

### EJ Success ###

EJSuccess <- chi.table(cases$Multinationals,cases$EJ.served.or.not)
EJSuccess_plot <- chi.plot(EJSuccess, cases$Multinationals)

### Rural vs Urban ###
rural <- chi.table(cases$Multinationals,cases$Type.of.Population)
rural_plot <- chi.plot(rural, cases$Multinationals)

### Actors mobilized ###

Actors <- chi.table(mobilizinggroups$Multinationals,mobilizinggroups)
Actors_plot <- chi.plot(Actors,mobilizinggroups$Multinationals)

### Forms of Mobilisation ###

Forms <- chi.table(mobilizingforms$Multinationals,mobilizingforms)
Forms_plot <- chi.plot(Forms,mobilizingforms$Multinationals)

### Env Impacts ###

env <- chi.table(envimpacts$Multinationals,envimpacts)
env_plot <- chi.plot(env,envimpacts$Multinationals)

### Hlt Impacts ###

hlt <- chi.table(hltimpacts$Multinationals,hltimpacts)
hlt_plot <- chi.plot(hlt,hltimpacts$Multinationals)

### Sec Impacts ###

sec <- chi.table(secimpacts$Multinationals,secimpacts)
sec_plot <- chi.plot(sec,secimpacts$Multinationals)

### Plots ###
Status_plot
Actors_plot
Forms_plot
Category_plot
Intensity_plot
EJSuccess_plot
rural_plot
env_plot
hlt_plot
sec_plot
Outcomes_plot
Types_plot
Commodities_plot

save.plot(Status_plot)
save.plot(Actors_plot)
save.plot(Forms_plot)
save.plot(Category_plot)
save.plot(Intensity_plot)
save.plot(EJSuccess_plot)
save.plot(rural_plot)
save.plot(env_plot)
save.plot(hlt_plot)
save.plot(sec_plot)
save.plot(Outcomes_plot)
save.plot(Types_plot)
save.plot(Commodities_plot)

#Save chitables
write.csv(Status,"Status_table.csv")
write.csv(Actors,"Actors_table.csv")
write.csv(Forms,"Forms_table.csv")
write.csv(Category,"Category_table.csv")
write.csv(Intensity,"Intensity_table.csv")
write.csv(EJSuccess,"EJSuccess_table.csv")
write.csv(rural,"rural_table.csv")
write.csv(env,"env_table.csv")
write.csv(hlt,"hlt_table.csv")
write.csv(sec,"sec_table.csv")
write.csv(Outcomes,"Outcomes_table.csv")
write.csv(Types,"Types_table.csv")
write.csv(Commodities,"Commodities_table.csv")

#### 25. Calculate whether env, hlt, and sec impacts are more common in case with foreign companies ####
#sum the total number of visible env, hlt, and sec impacts grouping by Multinationals
for(i in 1:nrow(cases)){
  cases$envimpacts.count[i] <- length(which(envimpacts[i,] %in% "V"))
  cases$hltimpacts.count[i] <- length(which(hltimpacts[i,] %in% "V"))
  cases$secimpacts.count[i] <- length(which(secimpacts[i,] %in% "V"))
}

env_count <- aggregate(cases$envimpacts.count, by=list(Category=cases$Multinationals), FUN=sum)
hlt_count <- aggregate(cases$hltimpacts.count, by=list(Category=cases$Multinationals), FUN=sum)
sec_count <- aggregate(cases$secimpacts.count, by=list(Category=cases$Multinationals), FUN=sum)

#Create data frame with total number of visible impacts with factors Multinationals and Type of impacts
No <- c(env_count[1,2],hlt_count[1,2],sec_count[1,2])
Yes <- c(env_count[2,2],hlt_count[2,2],sec_count[2,2])
Type <- c("Environmental","Health","Socioeconomic")
Impactschi <- cbind.data.frame(Type, Yes, No)
rm(env_count,hlt_count,sec_count,Yes,No,Type)

#Compute chi squared table
Impactschi$Ratio <- as.numeric(rep(NA,length.out = nrow(Impactschi)))
Residuals <-NULL
for(i in 1:nrow(Impactschi)){
  Impactschi$Ratio[i] <- Impactschi[i,2]/sum(Impactschi[i,2:3])
  chisq <- chisq.test(x = Impactschi[i,2:3], p = table(cases$Multinationals)[c(2,1)]/nrow(cases))
  Impactschi$p.value[i] <- chisq$p.value
  Res <- chisq$residuals
  Residuals <- rbind(Residuals,Res)
  Impactschi$significant[i] <- if(Impactschi$p.value[i] > 0.05) {"no"} else {"yes"}
}
Impactschi <- cbind(Impactschi,Residuals)
rm(Residuals,Res,chisq)

#Prepare data for ploting
Impactschi_N <- Impactschi[,-2]
Impactschi_M <- Impactschi[,-3]
colnames(Impactschi_M)[2] <- "data"
colnames(Impactschi_N)[2] <- "data"
Impactschi <- rbind(Impactschi_M,Impactschi_N)
rm(Impactschi_M,Impactschi_N)
Impactschi$Multinationals <- c("yes","yes","yes","no","no","no")
Impactschi$Type <- as.factor(Impactschi$Type)
levels(Impactschi$Type) <- c("Environmental (n = 11109)","Health (n = 3384)","Socioeconomic (n = 9069)")

#Plot
Impactschi_plot<-ggplot(Impactschi, aes(x = Type, y = data, fill=Multinationals)) +
  geom_col(position="fill",) +
  labs(x = NULL, y = NULL)  +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Foreign companies",values=c("#1e818b","#e2c537","#F6EB13","#7A752F","#734C20","#111111","#9B989B","#BE85BA","#3BB449","#ED2224","#AF3456"))+
  coord_flip()+
  geom_hline(yintercept = 1-(table(cases$Multinationals)[1]/nrow(cases)), size=0.5)

Impactschi_plot
save.plot(Impactschi_plot)
write.csv(Impactschi,"Impactschi_table.csv")

#### 26. How many conflicts are companies with over n conflicts each involved?  #### 

cases$count <- rep("no",nrow(cases))
Threshold <- 7
for(i in 1:nrow(cases)){
  if(any(Summary$Num_cases[grep(paste0("\\b",cases$Conflict.Id[i],"\\b"),Summary$conflictIDs)] > Threshold)){
    cases$count[i] <- "yes"
  }
} 
table(cases$count)
#Calculate for 0, 1, 4, 7 and 30 and plot in Adobe Illustrator
rm(Threshold)

#### 27. Extract top 96 superconflictive companies and assign UNGC level ####

Superconflictive <- Summary[1:96,]
Superconflictive$UNGC <- rep("",nrow(Superconflictive))
write.csv(Superconflictive, "Superconflictive.csv")
#Assign UNGC level manually in spreadsheet
Superconflictive <- read.csv("Superconflictive.csv")

table(Superconflictive$UNGC)
table(Superconflictive$Region)
table(Superconflictive$Category)
table(Superconflictive$Conflicts.overseas.)

#Returns the number of companies in each category in each UNGC level
X <- Superconflictive %>% group_by(Category,UNGC) %>% summarize(count=n())
X

#### 28. Where are the conflicts with companies from medium and low income countries? ####

MLISummary <- Summary[which(Summary$Region %in% c("Middle Income","Low Income")),]
table(MLISummary$Region)

MLISankey<-NULL
for(i in 1:nrow(MLISummary)){
  a <- as.numeric(unlist(strsplit(MLISummary$conflictIDs[i]," ")))
  Countries <- as.character(cases$Country[which(cases$Conflict.Id %in% a)])
  Companies <- rep(MLISummary$Country[i],length(Countries))
  d <- cbind.data.frame(Countries,Companies)
  MLISankey <- rbind(MLISankey,d)
}
rm(a,d, Countries, Companies)

DELETE <- NULL
for(i in 1:nrow(MLISankey)){
  if(Regions$World..Bank[which(Regions$Country %in% MLISankey[i,1])] %in% c("High Income","China"))
  {DELETE <- c(DELETE,i)}
}
MLISankey <- MLISankey[-DELETE,]

X<-0
for(i in 1:nrow(MLISankey)){
  if(MLISankey[i,1] == MLISankey[i,2]){X <- X+1}
}
X/nrow(MLISankey)
#87% of conflicts of companies from Middle and Low income countries are domestic
