# This script aims to investigate whether plants are consistent in their
# flowering phenology and whether there is population differentiation in
# consistency in FFDs between pop'ns.

# Uploading files
# setwd("C:\\Users\\lenglish\\Dropbox\\INBphenology") # for PC
setwd("~/Dropbox/INBphenology/consistencyReport/originalData") # Mac, dropbox data
library(echPhenology)
# Phen data all years 
p.05.15 <- read.csv("exPt1Phenology.csv")
p.10 <- read.csv("2010.CG1.Phenology.txt")
ped <- read.csv("96979899qGenPedigreeLE.csv")
rowpos <- read.csv("1996rowPosData.csv")

###### Editing dataframe, adding in data from 2010 ####

p.10$note <- NULL
names(p.10)[names(p.10) == "Year"] <- "year"
names(p.10)[names(p.10) == "startDateEarly"] <- "startDtEarly"
names(p.10)[names(p.10) == "startDateLate"] <- "startDtLate"
names(p.10)[names(p.10) == "endDateEarly"] <- "endDtEarly"
names(p.10)[names(p.10) == "endDateLate"] <- "endDtLate"
p.05.15$phenNote <- NULL
p.05.15$ttColor <- NULL
p.05.15$row <- NULL
p.05.15$pos <- NULL
p.05.15$cgHdId <- NULL

p.05.15$year <- substring(p.05.15$startDtEarly, 1, 4)
yy <-aggregate(startDtEarly~cgPlaId+year, data=p.05.15,length)
names(yy)[3] <- c("headCt")
p.05.15 <- merge(p.05.15, yy, by=c("cgPlaId", "year"))
p.05.15$phenYear <- NULL

# Formating dates to fit with the phenology package
p.05.15 <- updateDT(p.05.15, "startDtEarly", "sd.E")
p.05.15 <- updateDT(p.05.15, "startDtLate", "sd.L")
p.05.15 <- updateDT(p.05.15, "endDtEarly", "ed.E")
p.05.15 <- updateDT(p.05.15, "endDtLate", "ed.L")

p.10 <- updateDT(p.10, "startDtEarly", "sd.E")
p.10 <- updateDT(p.10, "startDtLate", "sd.L")
p.10 <- updateDT(p.10, "endDtEarly", "ed.E")
p.10 <- updateDT(p.10, "endDtLate", "ed.L")
p.10$startDtEarly <- NULL
p.10$startDtLate <- NULL
p.10$endDtEarly <- NULL
p.10$endDtLate <- NULL

# Aggregating phenolody data by plant instead of head and adding 2010 data 
class(p.05.15$sd.E)
st <- aggregate(p.05.15[ , c("sd.E", "sd.L")], 
                by = list(year = p.05.15$year, headCt =p.05.15$headCt, cgPlaId = p.05.15$cgPlaId), 
                FUN = min)
ed <- aggregate(p.05.15[ , c("ed.E", "ed.L")], 
                by = list(year = p.05.15$year, headCt = p.05.15$headCt, cgPlaId = p.05.15$cgPlaId), 
                FUN = max)
p.05.15.plant <- merge(st, ed)
## Phenology for all plants in CG1 ##
p.10 <- p.10[,c(2,3,1,4,5,6,7)]
pAllYrs <- rbind(p.05.15.plant, p.10[ ,names(p.05.15.plant)]) 

# removing records with unknown start date or bad start date
dim(pAllYrs[pAllYrs$sd.E < "1940-01-01", ])
dim(pAllYrs[pAllYrs$ed.L < "1940-01-01", ])

pAllYrs <- pAllYrs[pAllYrs$sd.E > "1940-01-01", ]
pAllYrs <- pAllYrs[pAllYrs$ed.L > "1940-01-01", ]

# Subsetting 1996 garden (try subsetting more gardens in the future?)
x1996 <- pAllYrs[pAllYrs$cgPlaId < 647, ] # 1602 records 

## adding column with Day of Year number
x1996$sn.E <- strptime(x1996$sd.E, format="%Y-%m-%d")$yday + 1
x1996$sn.L <- strptime(x1996$sd.L, format ="%Y-%m-%d")$yday + 1
x1996$en.E <- strptime(x1996$ed.E, format="%Y-%m-%d")$yday + 1
x1996$en.L <- strptime(x1996$ed.L, format ="%Y-%m-%d")$yday + 1

x1996 <- x1996[!x1996$sn.E == "1", ]

## adding column with date/month of FFD
x1996$sdm.E <- substr(x1996$sd.E, 6, 10)
x1996$sdm.L <- substr(x1996$sd.L, 6, 10)
x1996$edm.E <- substr(x1996$ed.E, 6, 10)
x1996$edm.L <- substr(x1996$ed.L, 6, 10)

# Adding a column with "phenCt" -- number of times an individual has flowered over the study period
nf <- aggregate(year ~ cgPlaId, x1996, function(x) length(unique(x))) 
names(nf)[names(nf)=="year"] <- "phenCt"
x1996 <- merge(x1996, nf, by.x="cgPlaId", by.y="cgPlaId")

# just getting phenCt
library(plyr)
ped <- read.csv("96979899qGenPedigreeLE.csv")
names(ped)[names(ped)=="cgplaid"] <- "cgPlaId"
ped96 <- ped[ped$cgPlaId < 647, ]
phenCts <- join(ped96, nf, type = "left")
phenCts <- phenCts[phenCts$site != "", ]
phenCts <- phenCts[phenCts$site != "Unknown", ]
phenCts <- phenCts[phenCts$site != "tower", ] # Only 1 plant flowered
newDF <- subset(phenCts, is.na(phenCts$phenCt)) # all bad plants - no phenCt
badPlas <- newDF[,"cgPlaId"]
phenCts <- phenCts[!(phenCts$cgPlaId %in% badPlas), ]
phenCts$site <- phenCts$site[ , drop=TRUE]
levels(phenCts$site)

x2015 <- x1996[x1996$year=="2015", ] #136 plants flowering that year
x2014 <- x1996[x1996$year=="2014", ] #93
x2013 <- x1996[x1996$year=="2013", ] #110 
x2012 <- x1996[x1996$year=="2012", ] #51 ""
x2011 <- x1996[x1996$year=="2011", ] #201 ""
x2010 <- x1996[x1996$year=="2010", ] #118 ""
x2009 <- x1996[x1996$year=="2009", ] #169 ""
x2008 <- x1996[x1996$year=="2008", ] #159 ""
x2007 <- x1996[x1996$year=="2007", ] #217 ""
x2006 <- x1996[x1996$year=="2006", ] #198 ""
x2005 <- x1996[x1996$year=="2005", ] #136 ""

# 3 ways of ranking FFD: DAF, DAM, and DOY ####
x1996 <- x1996[order(x1996[, "year"]), ]

## DAF ("Day after flowering") First day a plant flowers gets a rank of 1. Does take into account
## gaps in days in which no plant flowers. Example: if plant A starts flowering on day 1 and then 
## no plants flower for 3 days, plant B would start flowering on day 4 (not day 2)

min(x2005[,"sn.E"]) 
max(x2005[,"sn.E"])
or <- data.frame (x = 1:25, y = 176:200)
x2005 <- merge(x2005, or, by.x = "sn.E", by.y = "y")
names(x2005)[names(x2005)=="x"] <- "rank"

min(x2006[,"sn.E"]) 
max(x2006[,"sn.E"])
or <- data.frame (x = 1:31, y = 176:206)
x2006 <- merge(x2006, or, by.x = "sn.E", by.y = "y")
names(x2006)[names(x2006)=="x"] <- "rank"

min(x2007[,"sn.E"]) 
max(x2007[,"sn.E"])
or <- data.frame (x = 1:20, y = 171:190)
x2007 <- merge(x2007, or, by.x = "sn.E", by.y = "y")
names(x2007)[names(x2007)=="x"] <- "rank"

min(x2008[,"sn.E"]) 
max(x2008[,"sn.E"])
or <- data.frame (x = 1:37, y = 188:224)
x2008 <- merge(x2008, or, by.x = "sn.E", by.y = "y")
names(x2008)[names(x2008)=="x"] <- "rank"

min(x2009[,"sn.E"]) 
max(x2009[,"sn.E"])
or <- data.frame (x = 1:20, y = 179:198)
x2009 <- merge(x2009, or, by.x = "sn.E", by.y = "y")
names(x2009)[names(x2009)=="x"] <- "rank"

min(x2010[,"sn.E"]) 
max(x2010[,"sn.E"])
or <- data.frame (x = 1:25, y = 172:196)
x2010 <- merge(x2010, or, by.x = "sn.E", by.y = "y")
names(x2010)[names(x2010)=="x"] <- "rank"

min(x2011[,"sn.E"]) 
max(x2011[,"sn.E"])
or <- data.frame (x = 1:22, y = 186:207)
x2011 <- merge(x2011, or, by.x = "sn.E", by.y = "y")
names(x2011)[names(x2011)=="x"] <- "rank"

min(x2012[,"sn.E"]) 
max(x2012[,"sn.E"])
or <- data.frame (x = 1:24, y = 171:194)
x2012 <- merge(x2012, or, by.x = "sn.E", by.y = "y")
names(x2012)[names(x2012)=="x"] <- "rank"

min(x2013[,"sn.E"]) 
max(x2013[,"sn.E"])
or <- data.frame (x = 1:25, y = 189:213)
x2013 <- merge(x2013, or, by.x = "sn.E", by.y = "y")
names(x2013)[names(x2013)=="x"] <- "rank"

min(x2014[,"sn.E"]) 
max(x2014[,"sn.E"])
or <- data.frame (x = 1:20, y = 184:203)
x2014 <- merge(x2014, or, by.x = "sn.E", by.y = "y")
names(x2014)[names(x2014)=="x"] <- "rank"

min(x2015[,"sn.E"]) 
max(x2015[,"sn.E"])
or <- data.frame (x = 1:27, y = 184:210) 
x2015 <- merge(x2015, or, by.x = "sn.E", by.y = "y")
names(x2015)[names(x2015)=="x"] <- "rank"

x1996 <- do.call(rbind, list(x2005, x2006, x2007, x2008, x2009, x2010, x2011, x2012, x2013, x2014, x2015))

# Function for plotting plants on FFD graph
dd <- function(pla, color = "black") {
  lines(x1996[x1996$cgPlaId %in% pla, "year" ], x1996[x1996$cgPlaId %in% pla, "rank"], col = color, lwd = 2, type = "b", pch = 19)
}

plot(rank ~ year, x1996, xlab = "Year", ylab = "Day after First", pch = 1)
plot(x1996$year, jitter(x1996$rank), xlab = "Year", ylab = "Day after First", pch = 1)

## DAM ("Day after median") Ranks plants around central tendency for that given year (median first
## day of flowering = 0, three days before that = -3)


median(x2005$sn.E) # 185
min(x2005[,"sn.E"]) 
max(x2005[,"sn.E"])
rr <- data.frame (x = -9:15, y = 176:200)
x2005 <- merge(x2005, rr, by.x = "sn.E", by.y = "y")
names(x2005)[names(x2005)=="x"] <- "medRank"

median(x2006$sn.E) # 185 
min(x2006[,"sn.E"]) 
max(x2006[,"sn.E"])
rr <- data.frame (x = -9:21, y = 176:206)
x2006 <- merge(x2006, rr, by.x = "sn.E", by.y = "y")
names(x2006)[names(x2006)=="x"] <- "medRank"

median(x2007$sn.E) # 177
min(x2007[,"sn.E"]) 
max(x2007[,"sn.E"])
rr <- data.frame (x = -6:13, y = 171:190)
x2007 <- merge(x2007, rr, by.x = "sn.E", by.y = "y")
names(x2007)[names(x2007)=="x"] <- "medRank"

median(x2008$sn.E) #199
min(x2008[,"sn.E"]) 
max(x2008[,"sn.E"])
rr <- data.frame (x = -11:25, y = 188:224)
x2008 <- merge(x2008, rr, by.x = "sn.E", by.y = "y")
names(x2008)[names(x2008)=="x"] <- "medRank"

median(x2009$sn.E) # 185
min(x2009[,"sn.E"]) 
max(x2009[,"sn.E"])
rr <- data.frame (x = -6:13, y = 179:198)
x2009 <- merge(x2009, rr, by.x = "sn.E", by.y = "y")
names(x2009)[names(x2009)=="x"] <- "medRank"

median(x2010$sn.E) # 182
min(x2010[,"sn.E"]) 
max(x2010[,"sn.E"])
rr <- data.frame (x = -10:14, y = 172:196)
x2010 <- merge(x2010, rr, by.x = "sn.E", by.y = "y")
names(x2010)[names(x2010)=="x"] <- "medRank"

median(x2011$sn.E) # 195
min(x2011[,"sn.E"]) 
max(x2011[,"sn.E"])
rr <- data.frame (x = -9:12, y = 186:207)
x2011 <- merge(x2011, rr, by.x = "sn.E", by.y = "y")
names(x2011)[names(x2011)=="x"] <- "medRank"

median(x2012$sn.E) # 180
min(x2012[,"sn.E"]) 
max(x2012[,"sn.E"])
rr <- data.frame (x = -9:14, y = 171:194)
x2012 <- merge(x2012, rr, by.x = "sn.E", by.y = "y")
names(x2012)[names(x2012)=="x"] <- "medRank"

median(x2013$sn.E) # 200
min(x2013[,"sn.E"]) 
max(x2013[,"sn.E"])
rr <- data.frame (x = -11:13, y = 189:213)
x2013 <- merge(x2013, rr, by.x = "sn.E", by.y = "y")
names(x2013)[names(x2013)=="x"] <- "medRank"

median(x2014$sn.E) # 190
min(x2014[,"sn.E"]) 
max(x2014[,"sn.E"])
rr <- data.frame (x = -6:13, y = 184:203)
x2014 <- merge(x2014, rr, by.x = "sn.E", by.y = "y")
names(x2014)[names(x2014)=="x"] <- "medRank"

median(x2015$sn.E) # 198
min(x2015[,"sn.E"]) 
max(x2015[,"sn.E"])
rr <- data.frame (x = -14:12, y = 184:210)
x2015 <- merge(x2015, rr, by.x = "sn.E", by.y = "y")
names(x2015)[names(x2015)=="x"] <- "medRank"

x1996 <- do.call(rbind, list(x2005, x2006, x2007, x2008, x2009, x2010, x2011, x2012, x2013, x2014, x2015))

# plotting medRank
plot(medRank ~ year, x1996)
plot(x1996$year, jitter(x1996$medRank), xlab = "Year", ylab = "Day after Median",
     main = "", 
     abline(h=0, lty=3))

# Function to plot DAM (DM)
dm <- function(pla, color = "black") {
  lines(x1996[x1996$cgPlaId %in% pla, "year" ], x1996[x1996$cgPlaId %in% pla, "medRank"], col = color, lwd = 2, type = "b", pch = 19)
}

# plotting by median first date flowering

# 0 = median, 1 = 1 day after median, etc

## DOY ("Day of Year") Julian date plants began flowering. Especially useful for looking at 
## differences in timing between years
plot(x1996$year, jitter(x1996$sn.E), ylab = "Day of Year", xlab = "Year")
plot(x1996$year, jitter(x1996$sn.E), ylab = "Day of Year", xlab = "Year", yaxt="n",
     main="First flowering day distributions \n 1996 garden")
axis(2, at=c(170, 180, 190, 200, 210, 220), labels=c("Jun-19", "Jun-29", "July-9", "July-19", "July-29", "Aug-8"))
# Function for plotting DOY (DX)
dx <- function(pla, color = "black") {
  lines(x1996[x1996$cgPlaId %in% pla, "year" ], x1996[x1996$cgPlaId %in% pla, "sn.E"], col = color, lwd = 2, type = "b", pch = 19)
}

plot(en.E ~ year, x1996, ylab = "Day of Year", xlab = "Year")

#### Correlations between ranks between years ####

str(x1996)
plot(sn.E ~ en.E, x1996)
cor.test(x1996$sn.E, x1996$en.E, method = "spearman") # rho = .867
cor.test(x1996$sn.E, x1996$en.E, method = "pearson") # cor = .861
qqnorm(x1996$sn.E)
qqnorm(x1996$en.E)
x1996$dur <- x1996$en.L - x1996$sn.E
plot(sn.E ~ dur, x1996)

cgPlaId <- unique(x1996$cgPlaId)
new <- data.frame(cgPlaId, row.names = NULL)
x15 <- x1996[x1996$year=="2015",c("rank", "cgPlaId")]
# install.packages(plyr)
library(plyr)
test15 <- join(new, x15, type = "left")

x14 <- x1996[x1996$year=="2014",c("rank", "cgPlaId")]
x13 <- x1996[x1996$year=="2013",c("rank", "cgPlaId")]
x12 <- x1996[x1996$year=="2012",c("rank", "cgPlaId")]
x11 <- x1996[x1996$year=="2011",c("rank", "cgPlaId")]
x10 <- x1996[x1996$year=="2010",c("rank", "cgPlaId")]
x09 <- x1996[x1996$year=="2009",c("rank", "cgPlaId")]
x08 <- x1996[x1996$year=="2008",c("rank", "cgPlaId")]
x07 <- x1996[x1996$year=="2007",c("rank", "cgPlaId")]
x06 <- x1996[x1996$year=="2006",c("rank", "cgPlaId")]
x05 <- x1996[x1996$year=="2005",c("rank", "cgPlaId")]

test14 <- join(new, x14, type = "left")
test13 <- join(new, x13, type = "left")
test12 <- join(new, x12, type = "left")
test11 <- join(new, x11, type = "left")
test10 <- join(new, x10, type = "left")
test09 <- join(new, x09, type = "left")
test08 <- join(new, x08, type = "left")
test07 <- join(new, x07, type = "left")
test06 <- join(new, x06, type = "left")
test05 <- join(new, x05, type = "left")

tt <- do.call(cbind, list(test15, test14, test13, test12, test11, test10, test09, test08, test07, test06, test05))
tt[,3] <- NULL
tt[,4] <- NULL
tt[,5] <- NULL
tt[,6] <- NULL
tt[,7] <- NULL
tt[,8] <- NULL
tt[,9] <- NULL
tt[,10] <- NULL
tt[,11] <- NULL
tt[,12] <- NULL

# install.packages("data.table")
library(data.table)

names<- c("cgPlaId", "rank2015", "rank2014", "rank2013", "rank2012", "rank2011", "rank2010", "rank2009", 
          "rank2008","rank2007","rank2006","rank2005")
setnames(tt, names)

# Pairwise sample sizes

sum(complete.cases(tt$rank2015, tt$rank2014)) # 52
sum(complete.cases(tt$rank2015, tt$rank2013)) # 71
sum(complete.cases(tt$rank2015, tt$rank2012)) # 23
sum(complete.cases(tt$rank2015, tt$rank2011)) # 95
sum(complete.cases(tt$rank2015, tt$rank2010)) # 56
sum(complete.cases(tt$rank2015, tt$rank2009)) # 83
sum(complete.cases(tt$rank2015, tt$rank2008)) # 80
sum(complete.cases(tt$rank2015, tt$rank2007)) # 91
sum(complete.cases(tt$rank2015, tt$rank2006)) # 84
sum(complete.cases(tt$rank2015, tt$rank2005)) # 57

sum(complete.cases(tt$rank2014, tt$rank2013)) # 34
sum(complete.cases(tt$rank2014, tt$rank2012)) # 25
sum(complete.cases(tt$rank2014, tt$rank2011)) # 58
sum(complete.cases(tt$rank2014, tt$rank2010)) # 41
sum(complete.cases(tt$rank2014, tt$rank2009)) # 50
sum(complete.cases(tt$rank2014, tt$rank2008)) # 51
sum(complete.cases(tt$rank2014, tt$rank2007)) # 71
sum(complete.cases(tt$rank2014, tt$rank2006)) # 46
sum(complete.cases(tt$rank2014, tt$rank2005)) # 38

sum(complete.cases(tt$rank2013, tt$rank2012)) # 11
sum(complete.cases(tt$rank2013, tt$rank2011)) # 78
sum(complete.cases(tt$rank2013, tt$rank2010)) # 43
sum(complete.cases(tt$rank2013, tt$rank2009)) # 67
sum(complete.cases(tt$rank2013, tt$rank2008)) # 67
sum(complete.cases(tt$rank2013, tt$rank2007)) # 72
sum(complete.cases(tt$rank2013, tt$rank2006)) # 67
sum(complete.cases(tt$rank2013, tt$rank2005)) # 40

sum(complete.cases(tt$rank2011, tt$rank2012)) # 29
sum(complete.cases(tt$rank2010, tt$rank2012)) # 24
sum(complete.cases(tt$rank2009, tt$rank2012)) # 26
sum(complete.cases(tt$rank2008, tt$rank2012)) # 25
sum(complete.cases(tt$rank2007, tt$rank2012)) # 33
sum(complete.cases(tt$rank2006, tt$rank2012)) # 22
sum(complete.cases(tt$rank2005, tt$rank2012)) # 21

sum(complete.cases(tt$rank2011, tt$rank2010)) # 83
sum(complete.cases(tt$rank2011, tt$rank2009)) # 116
sum(complete.cases(tt$rank2011, tt$rank2008)) # 113
sum(complete.cases(tt$rank2011, tt$rank2007)) # 127
sum(complete.cases(tt$rank2011, tt$rank2006)) # 110
sum(complete.cases(tt$rank2011, tt$rank2005)) # 80

sum(complete.cases(tt$rank2010, tt$rank2009)) # 77
sum(complete.cases(tt$rank2010, tt$rank2008)) # 76
sum(complete.cases(tt$rank2010, tt$rank2007)) # 81
sum(complete.cases(tt$rank2010, tt$rank2006)) # 73
sum(complete.cases(tt$rank2010, tt$rank2005)) # 57

sum(complete.cases(tt$rank2008, tt$rank2009)) # 117
sum(complete.cases(tt$rank2007, tt$rank2009)) # 112
sum(complete.cases(tt$rank2006, tt$rank2009)) # 100
sum(complete.cases(tt$rank2005, tt$rank2009)) # 71

sum(complete.cases(tt$rank2007, tt$rank2008)) # 103
sum(complete.cases(tt$rank2006, tt$rank2008)) # 102
sum(complete.cases(tt$rank2005, tt$rank2008)) # 65

sum(complete.cases(tt$rank2006, tt$rank2007)) # 138
sum(complete.cases(tt$rank2005, tt$rank2007)) # 99

sum(complete.cases(tt$rank2005, tt$rank2006)) # 89

n15 <- c(NA, 52, 71, 23, 95, 56, 83, 80, 91, 84, 57)
n14 <- c(52, NA, 34, 25, 58, 41, 50, 51, 71, 46, 38)
n13 <- c(71, 34, NA, 11, 78, 43, 67, 67, 72, 67, 40)
n12 <- c(23, 25, 11, NA, 29, 24, 26, 25, 33, 22, 21)
n11 <- c(95, 58, 78, 29, NA, 83, 116, 113, 127, 110, 80 )
n10 <- c(56, 41, 43, 24, 83, NA, 77, 76, 81, 73, 57)
n09 <- c(83, 50, 67, 26, 116, 77, NA, 117, 112, 100, 71)
n08 <- c(80, 51, 67, 25, 113, 76, 117, NA, 103, 102, 65)
n07 <- c(91, 71, 72, 33, 127, 81, 112, 103, NA, 138, 99)
n06 <- c(84, 46, 67, 22, 110, 73, 100, 102, 138, NA, 89)
n05 <- c(57, 38, 40, 21, 80, 57, 71, 65, 99, 89, NA)

nn <- data.frame(year = c("2015", "2014", "2013", "2012", "2011", "2010", "2009", "2008", "2007", "2006", "2005"),
                 x2015 = n15, x2014 = n14, x2013=n13, x2012 = n12, x2011=n11, x2010=n10,x2009=n09, x2008=n08, x2007=n07, x2006=n06, x2005=n05)
stargazer(nn)

table <- cor(tt, use = "pairwise.complete.obs", method = "spearman")
plot(tt$rank2008, tt$rank2013, xlab = "Rank in 2008", ylab = "Rank in 2013", pch = 19)
text(30, 20, "p = 0.003")
plot(tt$rank2008, tt$rank2007, xlab = "Rank in 2008", ylab = "Rank in 2007", pch = 19)
text(30, 17, "p <0.001")
plot(tt$rank2008, tt$rank2009, xlab = "Rank in 2008", ylab = "Rank in 2009", pch = 19)
text(30, 17, "p = 0.95")

stargazer(tt)

cors <- data.frame(table)
cors[,"cgPlaId"] <- NULL
cors <- cors[-1,] 
cors

cor.test(tt$rank2015, tt$rank2014, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2013, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2010, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2015, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2014, tt$rank2013, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2010, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2014, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2010, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2013, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2006, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2005, tt$rank2012, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2006, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2005, tt$rank2011, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2010, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2012, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2011, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2009, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2008, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2007, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2006, exact = FALSE, method = "spearman")
cor.test(tt$rank2006, tt$rank2005, exact = FALSE, method = "spearman")

cor.test(tt$rank2013, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2012, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2011, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2010, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2009, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2008, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2007, tt$rank2005, exact = FALSE, method = "spearman")
cor.test(tt$rank2006, tt$rank2005, exact = FALSE, method = "spearman")

# p - value coorelation table 

p15 <- c("NA", "0.097","0.015*", "0.171", "0.001*", "0.063", "0.270","0.025*","0.015*","0.024*","0.181")
p14 <- c("0.097", "NA", "0.244", "0.103", "0.007*", "0.045*", "0.016*","0.886","0.016*","0.999","0.274")
p13 <- c("0.015*","0.244", "NA", "0.521", "0.004*", "0.172", "0.141","0.003*","0.010*","0.092","0.153")
p12 <- c("0.171","0.103", "0.521", "NA", "0.258" , "0.362","0.648","0.185","0.568","0.101", "0.167")
p11 <- c("0.001*", "0.007*", "0.004*", "0.258", "NA", "0.732", "0.086", "0.009*",  "0.001*", "0.009*",  "0.994")
p10 <- c("0.063", "0.045*", "0.172","0.362", "0.732", "NA", "0.006*", "0.004*", "0.034*", "0.015*", "0.215")
p09 <- c("0.270", "0.016*", "0.141", "0.648" , "0.086", "0.006*", "NA", "0.948", "<0.001*", "0.0467*", "0.305")
p08 <- c("0.025*", "0.886*", "0.003*", "0.185", "0.009*", "0.004*", "0.948", "NA", "0.001*", "0.001*","0.478")
p07 <- c("0.015*","0.016*", "0.010*", "0.568", "0.001*", "0.034*", "0.001*", "0.001*" ,"NA", "0.001*" ,"0.248")
p06 <- c("0.024*", "0.999", "0.092", "0.101", "0.009*", "0.015*" ,"0.0467*","0.001*","0.001*","NA","0.064")
p05 <- c("0.181", "0.274", "0.153", "0.167", "0.994", "0.215", "0.305", "0.478", "0.248", "0.064", "NA")

p15 <- c(NA, "0.097","0.015*", "0.171", "0.001*", "0.063", "0.270","0.025*","0.015*","0.024*","0.181")
p14 <- c(NA, NA, "0.244", "0.103", "0.007*", "0.045*", "0.016*","0.886","0.016*","0.999","0.274")
p13 <- c(NA,NA, NA, "0.521", "0.004*", "0.172", "0.141","0.003*","0.010*","0.092","0.153")
p12 <- c(NA, NA, NA, NA, "0.258" , "0.362","0.648","0.185","0.568","0.101", "0.167")
p11 <- c(NA, NA, NA, NA, NA, "0.732", "0.086", "0.009*",  "0.001*", "0.009*",  "0.994")
p10 <- c(NA, NA, NA,NA, NA, NA, "0.006*", "0.004*", "0.034*", "0.015*", "0.215")
p09 <- c(NA, NA, NA, NA , NA, NA, NA, "0.948", "<0.001*", "0.0467*", "0.305")
p08 <- c(NA, NA, NA, NA, NA, NA, NA, NA, "0.001*", "0.001*","0.478")
p07 <- c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, "0.001*" ,"0.248")
p06 <- c(NA, NA, NA, NA, NA, NA ,NA,NA,NA,NA,"0.064")
p05 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)




pp <- data.frame(year = c("2015", "2014", "2013", "2012", "2011", "2010", "2009", "2008", "2007", "2006", "2005"),
                 x2015 =p15, x2014=p14, x2013=p13, x2012 = p12, x2011=p11, x2010=p10,x2009=p09, x2008=p08, x2007=p07, x2006=p06, x2005=p05, stringsAsFactors = FALSE)

class(pp)
table(pp)
str(pp)
stargazer(pp, summary = FALSE) 

# rho correlation tables
r15 <- c("NA", "0.233", "0.288", "-0.296", "0.399","0.250","0.122","0.251","0.235", "0.246", "0.180")
r14 <- c("0.233", "NA", "0.205", "0.334", "0.351", "0.314","0.339","0.206","0.286","-0.001", "0.182")
r13 <- c("0.288", "0.205", "NA", "0.218", "0.323", "0.212", "0.182","0.360","0.303","0.209","-0.230")
r12 <- c("-0.296", "0.334", "0.218", "NA", "0.217" , "0.194","0.094","0.274","0.103","0.359", "0.313")
r11 <- c("0.399", "0.351", "0.323", "0.217", "NA", "0.038", "0.160", "0.246",  "0.414", "0.247",  "-0.001")
r10 <- c("0.250", "0.314", "0.212","0.194", "0.038", "NA", "0.309", "0.325", "0.236", "0.284", "0.167")
r09 <- c("0.122", "0.339", "0.182", "0.094" , "0.160", "0.309", "NA", "0.006", "0.322", "0.199", "0.124")
r08 <- c("0.251", "0.206", "0.360", "0.274", "0.246", "0.325", "0.006", "NA", "0.326", "0.427","0.090")
r07 <- c("0.235", "0.286", "0.303", "0.103", "0.414", "0.236", "0.322", "0.326" ,"NA", "0.332" ,"0.117")
r06 <- c("0.246", "-0.001", "0.209", "0.359", "0.247", "0.284" ,"0.199","0.427","0.332","NA","0.197")
r05 <- c("0.180", "0.182", "-0.230", "0.313", "-0.001", "0.167", "0.124", "0.090", "0.117", "0.197", "NA")

rr <- data.frame(year = c("2015", "2014", "2013", "2012", "2011", "2010", "2009", "2008", "2007", "2006", "2005"),
                 x2015=r15, x2014= r14, x2013=r13, x2012=r12, x2011=r11, x2010=r10, x2009=r09, x2008=r08, x2007=r07, x2006=r06, x2005=r05)

stargazer(rr, type = "text")



# Creating figure to visualize DAF, DOY, and DAM @ same time ####
x1996 <- x1996[order(x1996[, "year"]), ]

par(mfcol=c(3,3),
    oma = c(2,2,2,2) + 0.1,
    mar = c(4.1,4.1,1,1) + 0.1)


plot(rank ~ year, x1996, ylab = "Day after First", xlab = "")
dd("301", rainbow(15)[1])
dd("350", rainbow(15)[2])
dd("160", rainbow(15)[3])
dd("77", rainbow(15)[7])
dd("315", rainbow(15)[9])
dd("359", rainbow(15)[12])

plot(rank ~ year, x1996,  ylab = "Day after First", xlab = "")
dd("638", rainbow(15)[1])
dd("348", rainbow(15)[2])
dd("1", rainbow(15)[3])
dd("183", rainbow(15)[7])
dd("303", rainbow(15)[9])
dd("160", rainbow(15)[12])

plot(rank ~ year, x1996, xlab = "Year", ylab = "Day after First")
dd("40", rainbow(15)[1])
dd("52", rainbow(15)[2])
dd("61", rainbow(15)[3])
dd("66", rainbow(15)[7])
dd("81", rainbow(15)[9])
dd("89", rainbow(15)[12])

plot(sn.E ~ year, x1996, ylab = "Day of Year", xlab = "")
dx("301", rainbow(15)[1])
dx("350", rainbow(15)[2])
dx("160", rainbow(15)[3])
dx("77", rainbow(15)[7])
dx("315", rainbow(15)[9])
dx("359", rainbow(15)[12])

plot(sn.E ~ year, x1996, ylab = "Day of Year", xlab = "")
dx("638", rainbow(15)[1])
dx("348", rainbow(15)[2])
dx("1", rainbow(15)[3])
dx("183", rainbow(15)[7])
dx("303", rainbow(15)[9])
dx("160", rainbow(15)[12])

plot(sn.E ~ year, x1996, ylab = "Day of Year", xlab = "Year")
dx("40", rainbow(15)[1])
dx("52", rainbow(15)[2])
dx("61", rainbow(15)[3])
dx("66", rainbow(15)[7])
dx("81", rainbow(15)[9])
dx("89", rainbow(15)[12])

plot(medRank ~ year, x1996, xlab = "", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("301", rainbow(15)[1])
dm("350", rainbow(15)[2])
dm("160", rainbow(15)[3])
dm("77", rainbow(15)[7])
dm("315", rainbow(15)[9])
dm("359", rainbow(15)[12])

plot(medRank ~ year, x1996, xlab = "", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("638", rainbow(15)[1])
dm("348", rainbow(15)[2])
dm("1", rainbow(15)[3])
dm("183", rainbow(15)[7])
dm("303", rainbow(15)[9])
dm("160", rainbow(15)[12])

plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("40", rainbow(15)[1])
dm("52", rainbow(15)[2])
dm("61", rainbow(15)[3])
dm("66", rainbow(15)[7])
dm("81", rainbow(15)[9])
dm("89", rainbow(15)[12])

par(mfcol = c(1,1), 
    oma = c(2,2,2,2) + 0.1,
    mar = c(5.1,4.1,2.1,2.1) + 0.1)

# Creating a figure to visualize plants with phenCt = 8, 9, 10. Nine study years but 2012 is very low ###
x1996 <- x1996[order(x1996[, "year"]), ]
x10 <- x1996[x1996$phenCt==10, ]
x9 <- x1996[x1996$phenCt==9, ]
x8 <- x1996[x1996$phenCt==8, ]
length(unique(x8$cgPlaId)) # 24 plants have flowered 8 times
length(unique(x9$cgPlaId)) # 12 plants ""
length(unique(x10$cgPlaId)) # 4 plants
unique(x10$cgPlaId)
unique(x9$cgPlaId)
plot(medRank ~ year, x1996)
plot(x1996$year, jitter(x1996$medRank), xlab = "Year", ylab = "Day after Median",
     main = "", 
     abline(h=0, lty=3))

# plotting 4 plants that have flowered 10 times
dm <- function(pla, color = "black") {
  lines(x1996[x1996$cgPlaId %in% pla, "year" ], x1996[x1996$cgPlaId %in% pla, "medRank"], col = color, lwd = 2, type = "b", pch = 19)
}
plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
 dm("238", rainbow(15)[2])
 dm("610", rainbow(15)[9])
 dm("561", rainbow(15)[15])
 dm("278", rainbow(15)[5])
# 
plot(medRank ~ year, x10, xlab = "Year", ylab = "Day after Median", ylim=c(-10, 15))
abline(h = 0, lty = 3)
dm("238", rainbow(15)[2])
dm("610", rainbow(15)[9])
dm("561", rainbow(15)[13])
dm("278", rainbow(15)[8])

# plotting 12 plants that have flowered 9 times (too busy!)
# dm <- function(pla, color = "black") {
#   lines(x1996[x1996$cgPlaId %in% pla, "year" ], x1996[x1996$cgPlaId %in% pla, "medRank"], col = color, lwd = 2, type = "b", pch = 19)
# }
# plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
# abline(h = 0, lty = 3)
# dm("372", rainbow(15)[2])
# dm("323", rainbow(15)[9])
# dm("138", rainbow(15)[13])
# dm("401", rainbow(15)[8])
# dm("461", rainbow(15)[1])
# dm("359", rainbow(15)[5])
# dm("404", rainbow(15)[12])
# dm("30", rainbow(15)[15])
# dm("576", rainbow(15)[10])
# dm("225", rainbow(15)[11])
# dm("326", rainbow(15)[6])
# dm("152", rainbow(15)[3])
# 
# # 
# plot(medRank ~ year, x9, xlab = "Year", ylab = "Day after Median", ylim=c(-10, 15))
# abline(h = 0, lty = 3)
# dm("372", rainbow(15)[2])
# dm("323", rainbow(15)[9])
# dm("138", rainbow(15)[13])
# dm("401", rainbow(15)[8])
# dm("461", rainbow(15)[1])
# dm("359", rainbow(15)[5])
# dm("404", rainbow(15)[12])
# dm("30", rainbow(15)[15])
# dm("576", rainbow(15)[10])
# dm("225", rainbow(15)[11])
# dm("326", rainbow(15)[6])
# dm("152", rainbow(15)[3])

# 4 plots. Each one showing one of the 4 plants that flowered for all 10 years
par(mfcol=c(2, 2))
plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("238")
plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("610")
plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("561")
plot(medRank ~ year, x1996, xlab = "Year", ylab = "Day after Median")
abline(h = 0, lty = 3)
dm("278")
par(mfcol=c(1,1))

# Averaging DAM #### 

hist(nf$phenCt)
# Excluding plants that have only flowered once. 
mm <- x1996[x1996$phenCt != 1, ]
length(unique(mm$cgPlaId)) # n = 314

aggregate(medRank ~ cgPlaId, mm, length)
aggregate(medRank ~ cgPlaId, mm, range)
mns <- aggregate(medRank ~ cgPlaId, mm, mean)
names(mns)[names(mns)=="medRank"] <- "meanMedRank"
mm <- merge(mm, mns, by.x = "cgPlaId", by.y = "cgPlaId")
hist(mns$meanMedRank, 20, xlab = "Mean DAM Rank", main = "")
text(12, 30, "var = 8.52")
ts <- var(mns$meanMedRank)
mean(mns$meanMedRank)
mns <- mns[order(mns[, "meanMedRank"]),]
mm <- mm[order(mm[, "meanMedRank"]), ]

# DAM bootstrap analysis ####

mm <- mm[order(mm[,"year"]), ]
# 2005 = 1:129
# 2006 = 130:319
# 2007 = 320:529
# 2008 = 530:689
# 2009 = 690:857
# 2010 = 858:974
# 2011 = 975:1166 
# 2012 = 1167:1217  
# 2013 = 1218:1326
# 2014 = 1327:1419
# 2015 = 1420:1552

# An example of what we want the bootstrap funtion to do
bs.df <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
  (sample(mm[1:129, "medRank"], replace = TRUE)), 
  (sample(mm[130:319, "medRank"], replace = TRUE)),
  (sample(mm[320:529, "medRank"], replace = TRUE)),
  (sample(mm[530:689, "medRank"], replace = TRUE)),
  (sample(mm[690:857, "medRank"], replace = TRUE)),
  (sample(mm[858:974, "medRank"], replace = TRUE)),
  (sample(mm[975:1166, "medRank"], replace = TRUE)),
  (sample(mm[1167:1217, "medRank"], replace = TRUE)),
  (sample(mm[1218:1326, "medRank"], replace = TRUE)),
  (sample(mm[1327:1419, "medRank"], replace = TRUE)),
  (sample(mm[1420:1552, "medRank"], replace = TRUE))
))
bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
hist(bs.mns$bsMeanMedRank, 20)
var(bs.mns$bsMeanMedRank) # test statistic = 4.4

# trying to bootstrap

nboot <- 10000
bs.ts <- numeric(nboot)

for(i in 1:nboot){
  bs.df <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
    (sample(mm[1:129, "medRank"], replace = TRUE)), 
    (sample(mm[130:319, "medRank"], replace = TRUE)),
    (sample(mm[320:529, "medRank"], replace = TRUE)),
    (sample(mm[530:689, "medRank"], replace = TRUE)),
    (sample(mm[690:857, "medRank"], replace = TRUE)),
    (sample(mm[858:974, "medRank"], replace = TRUE)),
    (sample(mm[975:1166, "medRank"], replace = TRUE)),
    (sample(mm[1167:1217, "medRank"], replace = TRUE)),
    (sample(mm[1218:1326, "medRank"], replace = TRUE)),
    (sample(mm[1327:1419, "medRank"], replace = TRUE)),
    (sample(mm[1420:1552, "medRank"], replace = TRUE))
  ))
  bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
  names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
  bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
  bs.ts[i] <- var(bs.mns$bsMeanMedRank)
}
hist(bs.ts, 20, main = "", xlab = "Variance", xlim = c(3, 9))
mean(bs.ts)
abline(v=ts)

utpv <- length(bs.ts[bs.ts >= abs(ts)])/nboot
utpv                     
ltpv <- length(bs.ts[bs.ts <= abs(ts)])/nboot
ltpv                    
p <- 2*min(utpv, ltpv)   
p   # zero? less than 0.001

# bootstrap without outlier point, still highly significant 

# cgPlaId 350, flowers in 2008 and 2013
# mmTest <- mm[mm$cgPlaId != 350, ]
# 
# nboot <- 10000
# bs.ts <- numeric(nboot)
# 
# for(i in 1:nboot){
#   bs.df <- data.frame(cgPlaId = mmTest$cgPlaId, bsMedRank = c(
#     (sample(mmTest[1:129, "medRank"], replace = TRUE)), 
#     (sample(mmTest[130:319, "medRank"], replace = TRUE)),
#     (sample(mmTest[320:529, "medRank"], replace = TRUE)),
#     (sample(mmTest[530:688, "medRank"], replace = TRUE)),
#     (sample(mmTest[689:856, "medRank"], replace = TRUE)),
#     (sample(mmTest[857:973, "medRank"], replace = TRUE)),
#     (sample(mmTest[974:1165, "medRank"], replace = TRUE)),
#     (sample(mmTest[1166:1216, "medRank"], replace = TRUE)),
#     (sample(mmTest[1217:1324, "medRank"], replace = TRUE)),
#     (sample(mmTest[1325:1417, "medRank"], replace = TRUE)),
#     (sample(mmTest[1418:1550, "medRank"], replace = TRUE))
#   ))
#   bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
#   names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
#   bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
#   bs.ts[i] <- var(bs.mns$bsMeanMedRank)
# }
# hist(bs.ts, 20, main = "", xlab = "Variance", xlim = c(3, 9))
# mean(bs.ts)
# abline(v=ts)
# 
# utpv <- length(bs.ts[bs.ts >= abs(ts)])/nboot
# utpv                     
# ltpv <- length(bs.ts[bs.ts <= abs(ts)])/nboot
# ltpv                    
# p <- 2*min(utpv, ltpv)   
# p   # zero? less than 0.001

# # proportion of early vs. late plants (Not significant)
# plot(medRank ~ year, x1996)
# # early = more than 5 days before median, late = more than 5 days after
# ee <- length(mm[mm$medRank < -5, "medRank"])
# propE <- ee/1595 # 7.2% early
# ll <- length(mm[mm$medRank > 5, "medRank"])
# propL <- ll/1595 # 9.9% late
# 
# #proportion of plants who's average meanMedRank is above or below 5 and -5, respectively
# 
# nboot <- 10000
# probEb <- numeric(nboot)
# 
# for(i in 1:nboot){
#   bs.df3 <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
#     (sample(mm[1:129, "medRank"], replace = TRUE)), 
#     (sample(mm[130:319, "medRank"], replace = TRUE)),
#     (sample(mm[320:529, "medRank"], replace = TRUE)),
#     (sample(mm[530:689, "medRank"], replace = TRUE)),
#     (sample(mm[690:857, "medRank"], replace = TRUE)),
#     (sample(mm[858:974, "medRank"], replace = TRUE)),
#     (sample(mm[975:1166, "medRank"], replace = TRUE)),
#     (sample(mm[1167:1217, "medRank"], replace = TRUE)),
#     (sample(mm[1218:1326, "medRank"], replace = TRUE)),
#     (sample(mm[1327:1419, "medRank"], replace = TRUE)),
#     (sample(mm[1420:1552, "medRank"], replace = TRUE))
#   ))
#   eb <- length(bs.df3[bs.df3$bsMedRank < -5, "bsMedRank"])
#   probEb[i] <- eb/1552
# }
# hist(probEb, 20)

# nboot <- 10000
# probLb <- numeric(nboot)
# 
# for(i in 1:nboot){
#   bs.df4 <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
#     (sample(mm[1:129, "medRank"], replace = TRUE)), 
#     (sample(mm[130:319, "medRank"], replace = TRUE)),
#     (sample(mm[320:529, "medRank"], replace = TRUE)),
#     (sample(mm[530:689, "medRank"], replace = TRUE)),
#     (sample(mm[690:857, "medRank"], replace = TRUE)),
#     (sample(mm[858:974, "medRank"], replace = TRUE)),
#     (sample(mm[975:1166, "medRank"], replace = TRUE)),
#     (sample(mm[1167:1217, "medRank"], replace = TRUE)),
#     (sample(mm[1218:1326, "medRank"], replace = TRUE)),
#     (sample(mm[1327:1419, "medRank"], replace = TRUE)),
#     (sample(mm[1420:1552, "medRank"], replace = TRUE))
#   ))
#   lb <- length(bs.df4[bs.df4$bsMedRank > 5, "bsMedRank"])
#   probLb[i] <- lb/1552
# }
# hist(probLb, 20)

# extracting plants with average Dam of -1 to 1

middleDams <- subset(mm, meanMedRank >= -1 & meanMedRank <= 1) # n = 575 out of 1552 (little over a third)
middleDams <- middleDams[order(middleDams[,"year"]), ]
hist(middleDams$phenCt)
mean(middleDams$phenCt) # mean phenCt = 6.2 (slightly higher than all plants together)

# 2005 = 1:40
# 2006 = 41:100
# 2007 = 101:176
# 2008 = 177:237
# 2009 = 238:305
# 2010 = 306:355
# 2011 = 356:425 
# 2012 = 426:446  
# 2013 = 447:487
# 2014 = 488:527
# 2015 = 528:575

hist(middleDams$medRank, 20, xlab = "Day after Median scores", main = "")

#test run
middleMeans <- aggregate(medRank ~ cgPlaId, middleDams, mean) # 109 plants 
names(middleMeans)[names(middleMeans)=="medRank"] <- "meanMedRank"
hist(middleDams$medRank, 20, xlab = "Day after Median", 
main = "Day after median flowering scores for plants \n with an average DAM between -1 and 1")
hist(middleMeans$meanMedRank, 20, xlab = "Mean DAM Rank", main = "")
text(12, 30, "var = 0.47")
ts <- var(middleMeans$meanMedRank)

bs.df <- data.frame(cgPlaId = middleDams$cgPlaId, bsMedRank = c(
    (sample(middleDams[1:40, "medRank"], replace = TRUE)), 
    (sample(middleDams[41:100, "medRank"], replace = TRUE)),
    (sample(middleDams[101:176, "medRank"], replace = TRUE)),
    (sample(middleDams[177:237, "medRank"], replace = TRUE)),
    (sample(middleDams[238:305, "medRank"], replace = TRUE)),
    (sample(middleDams[306:355, "medRank"], replace = TRUE)),
    (sample(middleDams[356:425, "medRank"], replace = TRUE)),
    (sample(middleDams[426:446, "medRank"], replace = TRUE)),
    (sample(middleDams[447:487, "medRank"], replace = TRUE)),
    (sample(middleDams[488:527, "medRank"], replace = TRUE)),
    (sample(middleDams[528:575, "medRank"], replace = TRUE))
  ))
bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
hist(bs.mns$bsMeanMedRank, 20)
var(bs.mns$bsMeanMedRank) # example test = 2.17
  
# real deal
nboot <- 10000
bs.ts <- numeric(nboot)


for(i in 1:nboot){
  bs.df <- data.frame(cgPlaId = middleDams$cgPlaId, bsMedRank = c(
    (sample(middleDams[1:40, "medRank"], replace = TRUE)), 
    (sample(middleDams[41:100, "medRank"], replace = TRUE)),
    (sample(middleDams[101:176, "medRank"], replace = TRUE)),
    (sample(middleDams[177:237, "medRank"], replace = TRUE)),
    (sample(middleDams[238:305, "medRank"], replace = TRUE)),
    (sample(middleDams[306:355, "medRank"], replace = TRUE)),
    (sample(middleDams[356:425, "medRank"], replace = TRUE)),
    (sample(middleDams[426:446, "medRank"], replace = TRUE)),
    (sample(middleDams[447:487, "medRank"], replace = TRUE)),
    (sample(middleDams[488:527, "medRank"], replace = TRUE)),
    (sample(middleDams[528:575, "medRank"], replace = TRUE))
  ))
  bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
  names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
  bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
  bs.ts[i] <- var(bs.mns$bsMeanMedRank)
}
hist(bs.ts, 20, main = "", xlab = "Variance", xlim = c(0, 6))
abline(v=ts)
mean(bs.ts)

utpv <- length(bs.ts[bs.ts >= abs(ts)])/nboot
utpv                     
ltpv <- length(bs.ts[bs.ts <= abs(ts)])/nboot
ltpv                    
p <- 2*min(utpv, ltpv)   
p   # zero? less than 0.001

# bootstrap analysis for all plants with a DAM mean of less than -1 and greater than 1

outerDams.x <- subset(mm, meanMedRank > 1)
outerDams.y <- subset(mm, meanMedRank < -1)
outerDams <- rbind(outerDams.y, outerDams.x)
outerDams <- outerDams[order(outerDams[,"year"]), ]
hist(outerDams$medRank, 20)
hist(outerDams$phenCt)
mean(outerDams$phenCt) # mean phenCt = 5.63 (slightly higher than all plants together)

# 2005 = 1:89
# 2006 = 90:219
# 2007 = 220:353
# 2008 = 354:452
# 2009 = 453:552
# 2010 = 553:619
# 2011 = 620:741
# 2012 = 742:771  
# 2013 = 772:839
# 2014 = 840:892
# 2015 = 893:977

#test run
outerMeans <- aggregate(medRank ~ cgPlaId, outerDams, mean) # 205 plants 
names(outerMeans)[names(outerMeans)=="medRank"] <- "meanMedRank"
hist(outerDams$medRank, 20, xlab = "Day after Median", 
     main = "Day after median flowering scores for plants \n with an average DAM less than -1 and greater than 1")
hist(outerMeans$meanMedRank, 20, xlab = "Mean DAM Rank", main = "")
text(12, 30, "var = 12.53")
ts <- var(outerMeans$meanMedRank)

bs.df <- data.frame(cgPlaId = outerDams$cgPlaId, bsMedRank = c(
  (sample(outerDams[1:89, "medRank"], replace = TRUE)), 
  (sample(outerDams[90:219, "medRank"], replace = TRUE)),
  (sample(outerDams[220:353, "medRank"], replace = TRUE)),
  (sample(outerDams[354:452, "medRank"], replace = TRUE)),
  (sample(outerDams[453:552, "medRank"], replace = TRUE)),
  (sample(outerDams[553:619, "medRank"], replace = TRUE)),
  (sample(outerDams[620:741, "medRank"], replace = TRUE)),
  (sample(outerDams[742:771, "medRank"], replace = TRUE)),
  (sample(outerDams[772:839, "medRank"], replace = TRUE)),
  (sample(outerDams[840:892, "medRank"], replace = TRUE)),
  (sample(outerDams[893:977, "medRank"], replace = TRUE))
))
bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
hist(bs.mns$bsMeanMedRank, 20)
var(bs.mns$bsMeanMedRank) # example test = 5.744

# real deal
nboot <- 10000
bs.ts <- numeric(nboot)

for(i in 1:nboot){
  bs.df <- data.frame(cgPlaId = outerDams$cgPlaId, bsMedRank = c(
    (sample(outerDams[1:89, "medRank"], replace = TRUE)), 
    (sample(outerDams[90:219, "medRank"], replace = TRUE)),
    (sample(outerDams[220:353, "medRank"], replace = TRUE)),
    (sample(outerDams[354:452, "medRank"], replace = TRUE)),
    (sample(outerDams[453:552, "medRank"], replace = TRUE)),
    (sample(outerDams[553:619, "medRank"], replace = TRUE)),
    (sample(outerDams[620:741, "medRank"], replace = TRUE)),
    (sample(outerDams[742:771, "medRank"], replace = TRUE)),
    (sample(outerDams[772:839, "medRank"], replace = TRUE)),
    (sample(outerDams[840:892, "medRank"], replace = TRUE)),
    (sample(outerDams[893:977, "medRank"], replace = TRUE))
  ))
  bs.mns <- aggregate(bsMedRank ~ cgPlaId, bs.df, mean)
  names(bs.mns)[names(bs.mns)=="bsMedRank"] <- "bsMeanMedRank"
  bs.df <- merge(bs.df, bs.mns, by.x = "cgPlaId", by.y = "cgPlaId")
  bs.ts[i] <- var(bs.mns$bsMeanMedRank)
}
hist(bs.ts, 20, main = "", xlab = "Variance", xlim = c(2, 13))
abline(v=ts)
mean(bs.ts)

utpv <- length(bs.ts[bs.ts >= abs(ts)])/nboot
utpv                     
ltpv <- length(bs.ts[bs.ts <= abs(ts)])/nboot
ltpv                    
p <- 2*min(utpv, ltpv)   
p   # zero? less than 0.001

# Range bootstrap analysis ####

mm <- mm[order(mm[, "cgPlaId"]),] 
aggregate(medRank ~ cgPlaId, mm, range) # this works but only saves as two columns

# getting column with range
rx <- aggregate(medRank ~ cgPlaId, mm, max)
ry <- aggregate(medRank ~ cgPlaId, mm, min)
rz <- merge(rx, ry, by.x = "cgPlaId", by.y = "cgPlaId")
rz$range <- rz$medRank.x - rz$medRank.y # range of values
mm <- merge(mm, rz, by.x = "cgPlaId", by.y = "cgPlaId")
hist(mm$range, 20, xlab = "Per plant DAM range", main = "")
mean(mm$range)
ts2 <- var(mm$range)

# An example of what we want the bootstrap funtion to do
bs.df2 <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
  (sample(mm[1:129, "medRank"], replace = TRUE)), 
  (sample(mm[130:319, "medRank"], replace = TRUE)),
  (sample(mm[320:529, "medRank"], replace = TRUE)),
  (sample(mm[530:689, "medRank"], replace = TRUE)),
  (sample(mm[690:857, "medRank"], replace = TRUE)),
  (sample(mm[858:974, "medRank"], replace = TRUE)),
  (sample(mm[975:1166, "medRank"], replace = TRUE)),
  (sample(mm[1167:1217, "medRank"], replace = TRUE)),
  (sample(mm[1218:1326, "medRank"], replace = TRUE)),
  (sample(mm[1327:1419, "medRank"], replace = TRUE)),
  (sample(mm[1420:1552, "medRank"], replace = TRUE))
))
bx <- aggregate(bsMedRank ~ cgPlaId, bs.df2, max)
by <- aggregate(bsMedRank ~ cgPlaId, bs.df2, min)
bz <- merge(bx, by, by.x = "cgPlaId", by.y = "cgPlaId")
bz$bsRange <- bz$bsMedRank.x - bz$bsMedRank.y # range of values
bs.df2 <- merge(bs.df2, bz, by.x = "cgPlaId", by.y = "cgPlaId")
hist(bs.df2$bsRange, 20)
var(bs.df2$bsRange) # test statistic

nboot <- 10000
bs.ts <- numeric(nboot)

for(i in 1:nboot){
  bs.df2 <- data.frame(cgPlaId = mm$cgPlaId, bsMedRank = c(
    (sample(mm[1:129, "medRank"], replace = TRUE)), 
    (sample(mm[130:319, "medRank"], replace = TRUE)),
    (sample(mm[320:529, "medRank"], replace = TRUE)),
    (sample(mm[530:689, "medRank"], replace = TRUE)),
    (sample(mm[690:857, "medRank"], replace = TRUE)),
    (sample(mm[858:974, "medRank"], replace = TRUE)),
    (sample(mm[975:1166, "medRank"], replace = TRUE)),
    (sample(mm[1167:1217, "medRank"], replace = TRUE)),
    (sample(mm[1218:1326, "medRank"], replace = TRUE)),
    (sample(mm[1327:1419, "medRank"], replace = TRUE)),
    (sample(mm[1420:1552, "medRank"], replace = TRUE))
  ))
  bx <- aggregate(bsMedRank ~ cgPlaId, bs.df2, max)
  by <- aggregate(bsMedRank ~ cgPlaId, bs.df2, min)
  bz <- merge(bx, by, by.x = "cgPlaId", by.y = "cgPlaId")
  bz$bsRange <- bz$bsMedRank.x - bz$bsMedRank.y # range of values
  bs.df2 <- merge(bs.df2, bz, by.x = "cgPlaId", by.y = "cgPlaId")
  bs.ts[i] <- var(bs.df2$bsRange)
}
mean(bs.ts)
# p-value for range
utpv <- length(bs.ts[bs.ts >= abs(ts2)])/nboot
utpv                     
ltpv <- length(bs.ts[bs.ts <= abs(ts2)])/nboot
ltpv                    
p <- 2*min(utpv, ltpv)   
p   # 0.06, marginally significant   
hist(bs.ts, 20, xlab = "Variance", main = "")
abline(v=ts2)
text(30, 1200, "p=0.06")

# Adding in pedigree data ########
ped96 <- ped[ped$cgPlaId < 647, ]
names(ped96)[names(ped96)=="siteOfOriginPedigree"] <- "site"
names(ped96)[names(ped96)=="cgplaid"] <- "cgPlaId"
# fixing levels in pedigree data
levels(ped96$site)[levels(ped96$site) == "AA"] <- "aa"
levels(ped96$site)[levels(ped96$site) == "Eriley"] <- "eri"
levels(ped96$site)[levels(ped96$site) == "Lf"] <- "lf"
levels(ped96$site)[levels(ped96$site) == "NWLF"] <- "nwlf"
levels(ped96$site)[levels(ped96$site) == "SPP"] <- "spp"
levels(ped96$site)[levels(ped96$site) == "Stevens"] <- "sap"
levels(ped96$site)[levels(ped96$site) == "Nessman"] <- "ness"

# merging data
x1996 <- merge(x1996, ped96, by.x = "cgPlaId", by.y = "cgPlaId")

### get rid of empty levels ### 
# Only need to get rid of empty levels for site comparison analysis, can keep
# plants without pedigree data for consistency analysis
z1996 <- x1996[x1996$site != "", ]
z1996 <- z1996[z1996$site != "Unknown", ]
z1996 <- z1996[z1996$site != "Tower", ] # Only 1 plant flowered
z1996$site <- z1996$site[ , drop=TRUE]
levels(z1996$site)
### Using x1996 for consistency analysis and z1996 for site analysis ### 

# adding in row and position data, only need this for site comparison analysis 
z1996 <- merge(z1996, rowpos, by.x = "cgPlaId", by.y = "cgPlaId")

z2015 <- z1996[z1996$year=="2015", ]
z2014 <- z1996[z1996$year=="2014", ]
z2013 <- z1996[z1996$year=="2013", ]
z2012 <- z1996[z1996$year=="2012", ]
z2011 <- z1996[z1996$year=="2011", ]
z2010 <- z1996[z1996$year=="2010", ]
z2009 <- z1996[z1996$year=="2009", ]
z2008 <- z1996[z1996$year=="2008", ]
z2007 <- z1996[z1996$year=="2007", ]
z2006 <- z1996[z1996$year=="2006", ]
z2005 <- z1996[z1996$year=="2005", ]

# testing mean DAM values ####
# remove tower site - only 1 plant 
# remove plants that only flowered once

mz <- z1996[z1996$phenCt != 1, ]
aggregate(medRank ~ cgPlaId, mz, length)
means <- aggregate(medRank ~ cgPlaId, mz, mean) # aggregating cgPlaId by mean
names(means)[names(means)=="medRank"] <- "meanMedRank"
mz <- merge(mz, means, by.x = "cgPlaId", by.y = "cgPlaId")
plot(meanMedRank ~ site, mz) #meanMedRank by site
aggregate(meanMedRank ~ site, mz, mean)
str(mz)

# adding column for phenCt levels
mean(mz$phenCt)
mz$phenLevel <- "mid"
mz$phenLevel[mz$phenCt > 7] <- "high"
mz$phenLevel[mz$phenCt < 5] <- "low"
mz$phenLevel <- as.factor(mz$phenLevel)
levels(mz$phenLevel)
aggregate(cgPlaId ~ phenLevel, mz, length)

# visualizing distributions
hist(mz$meanMedRank)
hist(mz[mz$site == "aa", "meanMedRank"], 20)
hist(mz[mz$site == "eri", "meanMedRank"], 20)
hist(mz[mz$site == "lf", "meanMedRank"], 20)
hist(mz[mz$site == "ness", "meanMedRank"], 20)
hist(mz[mz$site == "nwlf", "meanMedRank"], 20)
hist(mz[mz$site == "spp", "meanMedRank"], 20)
hist(mz[mz$site == "sap", "meanMedRank"], 20)
str(mz)

m0 <- lm(meanMedRank ~ 1, mz)
m1 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row + pos + row:pos, mz) # no
m2 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row + pos, mz) # yes
m3 <- lm(meanMedRank ~ site + phenCt + row + pos, mz) # no
m4 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row, mz) # no
m5 <- lm(meanMedRank ~ site + phenCt + site:phenCt + pos, mz) # yes
m6 <- lm(meanMedRank ~ site + phenCt + site:phenCt, mz)
m7 <- lm(meanMedRank ~ site + phenCt + pos, mz)
m8 <- lm(meanMedRank ~ site + pos, mz)

str(mz)

anova(m1, m2) # Not significantly different, go with m2
anova(m2, m3)# Yes sig. diff. Interactive site/phenCt term matters. Start with m2
anova(m4, m2) # Yes sig. diff. Pos matters
anova(m5, m2) # Not sig. diff. Row doesn't matter. 
anova(m6, m5) # Pos matters
anova(m7, m5) # interactive term matters
anova(m0, m5)
anova(m7, m8)

hist(mz$phenCt)

summary (m2)

### here's our model m2  - meanMedRank ~ site + phenCt + site:phenCt + pos + row

# BUT
m8 <- lm(meanMedRank ~ row + pos + row:pos, mz)
m9 <- lm(meanMedRank ~ row + pos, mz)
m10 <- lm(meanMedRank ~ pos, mz)
anova(m8, m9) # no sig difference, interactive term doesn't matter
anova(m9, m10) # Yes sig difference, row matters when looking at just row and pos

par(mfcol = c(2, 2))
plot(m2)  # looks OK?
par(mfcol = c(1, 1))

plot(meanMedRank ~ phenCt, mz)
a2 <- mz[mz$phenCt==2, ]
a3 <- mz[mz$phenCt==3, ]
a4 <- mz[mz$phenCt==4, ]
a5 <- mz[mz$phenCt==5, ]
a6 <- mz[mz$phenCt==6, ]
a7 <- mz[mz$phenCt==7, ]
a8 <- mz[mz$phenCt==8, ]

plot(meanMedRank ~ site, a2)
plot(meanMedRank ~ site, a3)
plot(meanMedRank ~ site, a4)
plot(meanMedRank ~ site, a5)
plot(meanMedRank ~ site, a6)
plot(meanMedRank ~ site, a7)
plot(meanMedRank ~ site, a8)

coef(m2) # THIS IS CONFUSING?!
summary(m2)

# predicted values
str(mz)

aa <- mz[mz$site=="aa", ]
eri <- mz[mz$site=="eri", ]
lf <- mz[mz$site=="lf", ]
nwlf <- mz[mz$site=="nwlf", ]
ness <- mz[mz$site=="ness", ]
sap <- mz[mz$site=="sap", ]
spp <- mz[mz$site=="spp", ]

range(mz$phenCt)
pv <- 2:10
pv <- rep(pv, 7)

pv2 <- c(rep("ness", 9), rep("lf", 9), rep("sap", 9), 
         rep("spp", 9), rep("aa", 9), rep("nwlf", 9), 
         rep("eri", 9))
pv2
newd <- data.frame(phenCt= pv, site= pv2)
newd$pos <- 950
newd$row <- 26
str(newd)
str(mz)
median(mz$row)
newd$fit <- predict(m2, newd, type = "response") # Using m6 bc doesn't have position
str(newd)

range(mz$meanMedRank)

plot(1, type="n", xlab="phenCount", ylab="mean DAM rank", xlim=c(2,10), ylim=c(-9, 20), 
     main = "Mean DAM of plants by the number of times they flowered over \n the study period. Colors indicate origin site and lines are predicted \n values of mean DAM for each site at each phenology count.")
points(nwlf$phenCt, nwlf$meanMedRank, pch = 20, col = "pink")
points(eri$phenCt, jitter(eri$meanMedRank), pch = 20, col = "purple")
points(spp$phenCt, jitter(spp$meanMedRank), pch = 20, col = "red")
points(ness$phenCt, jitter(ness$meanMedRank), pch = 20, col = "green")
points(lf$phenCt, lf$meanMedRank, pch = 20, col = "yellow")
points(sap$phenCt, sap$meanMedRank, pch = 20, col = "gray")
points(aa$phenCt, aa$meanMedRank, pch = 20, col = "blue")

lines(newd$phenCt[newd$site == "aa"],
      newd$fit[newd$site == "aa"],
      col = "blue")
lines(newd$phenCt[newd$site == "eri"],
      newd$fit[newd$site == "eri"],
      col = "purple")
lines(newd$phenCt[newd$site == "ness"],
      newd$fit[newd$site == "ness"],
      col = "green")
lines(newd$phenCt[newd$site == "lf"],
      newd$fit[newd$site == "lf"],
      col = "yellow")
lines(newd$phenCt[newd$site == "nwlf"],
      newd$fit[newd$site == "nwlf"],
      col = "pink")
lines(newd$phenCt[newd$site == "spp"],
      newd$fit[newd$site == "spp"],
      col = "red")
lines(newd$phenCt[newd$site == "sap"],
      newd$fit[newd$site == "sap"],
      col = "gray")

legend2 <- c("aa", "eri", "lf", "ness","nwlf", "sap", "spp")
legend('topright', legend2, 
       lty=1, col=c('blue','purple','yellow','green','pink','gray','red'), bty='n', cex=.75)


eek <- mz[!duplicated(mz$cgPlaId),]
aggregate(cgPlaId~phenCt+site, eek, length)


highs <- subset(mz, phenLevel == "high")
mids <- subset(mz, phenLevel == "mid")
lows <- subset(mz, phenLevel == "low")
highPlot <- aggregate(medRank~site, highs, mean)
addRow <- c("ness","")
highPlot <- rbind(highPlot,addRow)
highPlot <- highPlot[order(highPlot$site),]
highPlot$order <- c(1:7)
midPlot <- aggregate(medRank~site, mids, mean)
midPlot$order <- c(1:7)
lowPlot <- aggregate(medRank~site, lows, mean)
lowPlot$order <- c(1:7)

#range(highPlot$medRank) # -1.26 to 0.588
#range(midPlot$medRank) # -1.23 to 1.65
#range(lowPlot$medRank) # -1.76 to 3.29

# set y-axis from -2 to 3.5 

plot(1, type="n", xlab="Site", ylab="mean DAM rank", xlim=c(1,7), xaxt="n", ylim=c(-2, 3.5), 
     main = "")
axis(1, at=1:7, labels=c("aa", "eri", "lf", "ness", "nwlf", "spp", "sap"))
points(highPlot$site, highPlot$medRank, pch = 19, col = "blue")
points(midPlot$site, midPlot$medRank, pch = 15, col = "gold")
points(lowPlot$site, lowPlot$medRank, pch = 17, col = "red")
lines(midPlot$site, midPlot$medRank, col = "gold")
lines(highPlot$site, highPlot$medRank, col = "blue")
lines(lowPlot$site, lowPlot$medRank, col = "red")
legend <- c("high(flowered 8-10 times)", "mid(flowered 5-7 times)", "low(flowered 2-4 times)")
legend('topright', legend, 
       lty=1, col=c('blue', 'gold',' red'), bty='n', cex=.75)



# now trying as with phenCt as level 
# 
# m0 <- lm(meanMedRank ~ 1, mz)
# mA <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row + pos + row:pos, mz) # no
# mB <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row + pos, mz) # yes
# mC <- lm(meanMedRank ~ site + phenLevel + row + pos, mz) # no
# mD <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row, mz) # no
# mE <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + pos, mz) # yes
# mF <- lm(meanMedRank ~ site + phenLevel + site:phenLevel, mz)
# mG <- lm(meanMedRank ~ site + phenLevel + pos, mz)
# 
# anova(mA, mB) # Not significantly different, go with m2
# anova(mB, mC) # Yes sig. diff. Interactive site/phenCt term matters. Start with m2
# anova(mD, mB) # Yes sig. diff. Pos matters
# anova(mE, mB) # Not sig. diff. Row doesn't matter. 
# anova(mF, mE) # Pos matters
# anova(mG, mE) # interactive term matters
# anova(m0, mE) # End up with same model

# testing ranges per site ####

rx <- aggregate(medRank ~ cgPlaId, mz, max)
ry <- aggregate(medRank ~ cgPlaId, mz, min)
rz <- merge(rx, ry, by.x = "cgPlaId", by.y = "cgPlaId")
rz$range <- rz$medRank.x - rz$medRank.y # range of values
mz <- merge(mz, rz, by.x = "cgPlaId", by.y = "cgPlaId")
hist(mz$range, 20, xlab = "Per plant DAM range", main = "")
str(mz)

plot(range ~ site, mz)
aggregate(range ~ site, mz, mean)

nm <- lm(range ~ 1, mz)
am1 <- lm(range ~ site + phenCt + site:phenCt + row + pos + row:pos, mz) # no
am2 <- lm(range ~ site + phenCt + site:phenCt + row + pos, mz) # yes
am3 <- lm(range ~ site + phenCt + row + pos, mz) # 
am4 <- lm(range ~ site + phenCt + site:phenCt + row, mz) # no
am5 <- lm(range ~ site + phenCt + site:phenCt + pos, mz) # yes
am6 <- lm(range ~ site + phenCt + site:phenCt, mz)
am7 <- lm(range ~ site + phenLevel + pos, mz)
am8 <- lm(range ~ site, mz)
am9 <- lm(range ~ pos, mz)

anova(am1, am2)
anova(am2, am3) # site, phenCt matters
anova(am2, am4) # row matters
anova(am2, am5) # pos matters
anova(am2, nm)
anova(am9, nm)

# model = range ~ site + phenCt + site:phenCt + row + pos

aa <- mz[mz$site=="aa", ]
eri <- mz[mz$site=="eri", ]
lf <- mz[mz$site=="lf", ]
nwlf <- mz[mz$site=="nwlf", ]
ness <- mz[mz$site=="ness", ]
sap <- mz[mz$site=="sap", ]
spp <- mz[mz$site=="spp", ]

range(mz$phenCt)
pv <- 2:10
pv <- rep(pv, 7)

pv2 <- c(rep("ness", 9), rep("lf", 9), rep("sap", 9), 
         rep("spp", 9), rep("aa", 9), rep("nwlf", 9), 
         rep("eri", 9))
pv2
median(mz$row)
newd <- data.frame(phenCt= pv, site= pv2)
newd$pos <- 950
newd$row <- 26
str(newd)
str(mz)
newd$fit <- predict(am2, newd, type = "response")
str(newd)

plot(mz$phenCt, mz$range, type = "n", ylim=c(0,25), xlim = c(2,12), xlab = "phenCount", ylab = "DAM range")
points(aa$phenCt, aa$range, pch = 20, col = "blue")
points(eri$phenCt, eri$range, pch = 20, col = "purple")
points(ness$phenCt, ness$range, pch = 20, col = "green")
points(lf$phenCt, lf$range, pch = 20, col = "yellow")
points(nwlf$phenCt, nwlf$range, pch = 20, col = "pink")
points(spp$phenCt, spp$range, pch = 20, col = "red")
points(sap$phenCt, sap$range, pch = 20, col = "gray")

lines(newd$phenCt[newd$site == "aa"],
      newd$fit[newd$site == "aa"],
      col = "blue")
lines(newd$phenCt[newd$site == "eri"],
      newd$fit[newd$site == "eri"],
      col = "green")
lines(newd$phenCt[newd$site == "ness"],
      newd$fit[newd$site == "ness"],
      col = "purple")
lines(newd$phenCt[newd$site == "lf"],
      newd$fit[newd$site == "lf"],
      col = "red")
lines(newd$phenCt[newd$site == "nwlf"],
      newd$fit[newd$site == "nwlf"],
      col = "gold")
lines(newd$phenCt[newd$site == "spp"],
      newd$fit[newd$site == "spp"],
      col = "pink")
lines(newd$phenCt[newd$site == "sap"],
      newd$fit[newd$site == "sap"],
      col = "gray")

legend2 <- c("aa", "eri", "lf", "ness","nwlf", "sap", "spp")
legend('topright', legend2, 
       lty=1, col=c('blue','purple','yellow','green','pink','gray','red'), bty='n', cex=.75)


par(mfcol = c(2,2))
plot(m2)
par(mfcol = c(1,1))

# graph for range
highs2 <- subset(mz, phenLevel == "high")
mids2 <- subset(mz, phenLevel == "mid")
lows2 <- subset(mz, phenLevel == "low")
highPlot2 <- aggregate(range~site, highs2, mean)
addRow <- c("ness","")
highPlot2 <- rbind(highPlot2,addRow)
highPlot2 <- highPlot2[order(highPlot2$site),]
highPlot2$order <- c(1:7)
midPlot2 <- aggregate(range~site, mids2, mean)
midPlot2$order <- c(1:7)
lowPlot2 <- aggregate(range~site, lows2, mean)
lowPlot2$order <- c(1:7)
plot(1, type="n", xlab="Site", ylab="range of DAM", xlim=c(1,7), xaxt="n", ylim=c(0, 18), 
     main = "Range of DAM for plants from each of the origin sites \n Plants are further subdivided by number of times they \n flowered over the study period")
axis(1, at=1:7, labels=c("aa", "eri", "lf", "ness", "nwlf", "spp", "sap"))
points(highPlot2$site, highPlot2$range, pch = 19, col = "blue")
points(midPlot2$site, midPlot2$range, pch = 15, col = "gold")
points(lowPlot2$site, lowPlot2$range, pch = 17, col = "red")
lines(midPlot2$site, midPlot2$range, col = "gold")
lines(highPlot2$site, highPlot2$range, col = "blue")
lines(lowPlot2$site, lowPlot2$range, col = "red")
legend <- c("high(flowered 8-10 times)", "mid(flowered 5-7 times)", "low(flowered 2-4 times)")
legend('topright', legend, 
       lty=1, col=c('blue', 'gold',' red'), bty='n', cex=.75)

range(mz$range)

# nm <- lm(range ~ 1, mm)
# amA <- lm(range ~ site + phenLevel + site:phenLevel + row + pos + row:pos, mm)
# amB <- lm(range ~ site + phenLevel + site:phenLevel + row + pos, mm)
# amC <- lm(range ~ site + phenLevel + row + pos, mm)
# amD <- lm(range ~ site + phenLevel + site:phenLevel + row, mm)
# amE <- lm(range ~ site + phenLevel + site:phenLevel + pos, mm)
# amF <- lm(range ~ site + phenLevel + site:phenLevel, mm)
# amG <- lm(range ~ site + phenLevel + pos, mm)
# amH <- lm(range ~ site, mm)
# amI <- lm(range ~ pos, mm)
# 
# 
# anova(amA, amB) # N.S.
# anova(amB, amC) # sig. phenLevel matters
# anova(amB, amD)
# anova(amB, amE)

# Analyzing whether site influence phenCt ####
aggregate(cgPlaId ~ phenCt + site, phenCts, length)
aggregate(phenCt ~ site, phenCts, length)
aggregate(site ~ phenCt, phenCts, length)
hist(phenCts$phenCt,breaks=rep(1:8,each=2)+c(-.4,.4), freq = TRUE, xlim = c(0,9))
aggregate(phenCt ~ site, phenCts, mean)
median(phenCts$phenCt)
mean(phenCts$phenCt)

null <- lm(phenCt ~ 1, phenCts)
altM <- lm(phenCt ~ site, phenCts)
altM <- lm(phenCt ~ site, phenCts)

anova(null, altM) # Not significant. Site not a good predictor of number of times
# a plant flowers.


hist(phenCts[phenCts$site=="sap", "phenCt"])
hist(phenCts[phenCts$site=="ness", "phenCt"], 10)
hist(phenCts[phenCts$site=="ness", "phenCt"], 10)


x2013 <- x1996[x1996$year=="2013", ]
x2012 <- x1996[x1996$year=="2012", ]
x2011 <- x1996[x1996$year=="2011", ]
x2010 <- x1996[x1996$year=="2010", ]
x2009 <- x1996[x1996$year=="2009", ]
x2008 <- x1996[x1996$year=="2008", ]
x2007 <- x1996[x1996$year=="2007", ]
x2006 <- x1996[x1996$year=="2006", ]
x2005 <- x1996[x1996$year=="2005", ]

# basic info tables ######

#phenCt distribution
hist(x1996$phenCt, xlab = "Phenology Count", breaks = 10, 
     main = "Number of times plants flowered over \n the study period (2005-2015)")
     
# number of plants flowering each year
n.values <- aggregate(sd.E ~ year, x1996, length)
names(n.values)[names(n.values) == "sd.E"] <- "n"
stargazer(n.values, summary = FALSE)
# number of plants from each site in the garden 
siteNums<- aggregate(expNm ~ site, ped96, length)
names(siteNums)[names(siteNums) == "site"] <- "OriginSite"
names(siteNums)[names(siteNums) == "expNm"] <- "n"
siteNums <- siteNums[-c(1,10), ] 
stargazer(siteNums, summary = FALSE)

# death rate
cc <- read.csv("http://echinaceaproject.org/data/cg1CoreData.csv")
c96 <- cc[cc$yrPlanted == 1996, ]

sum(c96$ld2015) #304
sum(c96$ld2014) #318
sum(c96$ld2013) #329
sum(c96$ld2012) #334
sum(c96$ld2011) #339
sum(c96$ld2010) #346
sum(c96$ld2009) #348
sum(c96$ld2008) #366
sum(c96$ld2007) #398
sum(c96$ld2006) #405
sum(c96$ld2005) #422
sum(c96$ld2004) #423
sum(c96$ld2003) #437
sum(c96$ld2002) #459
sum(c96$ld2001) #508
sum(c96$ld2000) #540
sum(c96$ld1999) #571
sum(c96$ld1998) #594
sum(c96$ld1997) #606
sum(c96$ld1996) #646

death <- data.frame(year=1996:2015, n=c(646,606,594,571,540,508,459,437,423,422,405,398,366,348,346,339,334,329,318,304))
plot(death,ylab="Number of living plants", ylim=c(200,700), pch = 20, xaxt = "n")
axis(1, at=c(1996, 2001, 2006, 2011, 2015))
