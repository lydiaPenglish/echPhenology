# This script aims to investigate the consistency in FFDs between pop'ns. 

# setwd("C:\\Users\\lenglish\\Dropbox\\INBphenology")
setwd("~/Dropbox/INBphenology") # for mac
library(echPhenology)
p13 <- read.csv(url("http://dl.dropboxusercontent.com/s/45qxicxev44se3y/2013.cg1.phenology.csv")) # 2013 data phenology data
p.05.12 <- read.csv("http://dl.dropboxusercontent.com/s/xogapzviqlf5ec5/allYearsNot2010.CG1.Phenology.txt", stringsAsFactors = FALSE) # 2005-2012 phenology data, no 2010
p.10 <- read.csv("http://dl.dropboxusercontent.com/s/v1dpnrt0wklfttg/2010.CG1.Phenology.txt") # 2010 data
ped <- read.csv("http://dl.dropboxusercontent.com/s/cqmeo58g298i2xx/96979899qGenPedigreeLE.csv")
rowpos <- read.csv("http://dl.dropboxusercontent.com/s/uygb7qmxib3w3bf/1996rowPosData.csv")

p13 <- read.csv("")

###### editing and putting together dataframes ####
### altering phenology dataframes to be able to merge them 
p.10$note <- NULL
p.05.12$note <- NULL
p13$phenNote <- NULL
names(p.10)[names(p.10) == "Year"] <- "year"
names(p13)[names(p13)=="phenYear"] <- "year"
names(p.05.12)[names(p.05.12)=="startDateEarly"] <- "startDtEarly"
names(p.05.12)[names(p.05.12)=="startDateLate"] <- "startDtLate"
names(p.05.12)[names(p.05.12)=="endDateEarly"] <- "endDtEarly"
names(p.05.12)[names(p.05.12)=="endDateLate"] <- "endDtLate"
names(p.05.12)[names(p.05.12)=="cgHeadId"] <- "cgHdId"

# 2010 data already has head count for each plant. I need to calculate head
# count for 2013 and 2005-2012 data and add a column to the dataframe.
### 2005-2012
yy <-aggregate(startDtEarly~cgPlaId+year, data=p.05.12,length)
names(yy)[3] <- c("headCt")
p.05.12 <- merge(p.05.12, yy, by=c("cgPlaId", "year"))
### 2013
qq <- aggregate(startDtEarly ~ cgPlaId + year, data = p13, length)
names(qq)[3] <- c("headCt")
p13 <- merge (qq, p13, by = c("cgPlaId", "year")) 

## Putting together 2005-2013 data (minus 2010)
ptemp <- rbind(p.05.12, p13[ , names(p.05.12)])
# Aggregating phenolody data by plant instead of head and adding 2010 data 
st <- aggregate(ptemp[ , c("startDtEarly", "startDtLate")], 
                by = list(year = ptemp$year, headCt =ptemp$headCt, cgPlaId = ptemp$cgPlaId), 
                FUN = min)
ed <- aggregate(ptemp[ , c("endDtEarly", "endDtLate")], 
                by = list(year = ptemp$year, headCt = ptemp$headCt, cgPlaId = ptemp$cgPlaId), 
                FUN = max)
ptemp.plant <- merge(st, ed)

p.plant <- rbind(ptemp.plant, p.10[, names(ptemp.plant)]) #Phenology for all plants in plot

# Formating dates to fit with the phenology package
p.plant <- updateDT(p.plant, "startDtEarly", "sd.E")
p.plant <- updateDT(p.plant, "startDtLate", "sd.L")
p.plant <- updateDT(p.plant, "endDtEarly", "ed.E")
p.plant <- updateDT(p.plant, "endDtLate", "ed.L")

# removing records with unknown start date
dim(p.plant[p.plant$sd.E < "1940-01-01", ])
dim(p.plant[p.plant$ed.L < "1940-01-01", ])

p.plant <- p.plant[p.plant$sd.E > "1940-01-01", ]
p.plant <- p.plant[p.plant$ed.L > "1940-01-01", ]

# reference pop
ALL1990s <- p.plant[p.plant$cgPlaId < 2500, ] 

# 1996
x1996 <- p.plant[p.plant$cgPlaId < 647, ] # phen data for all plants in the 1996 garden

## adding column to ss with Day of Year number
x1996$sn.E <- strptime(x1996$sd.E, format="%Y-%m-%d")$yday + 1
x1996$sn.L <- strptime(x1996$sd.L, format ="%Y-%m-%d")$yday + 1
x1996$en.E <- strptime(x1996$ed.E, format="%Y-%m-%d")$yday + 1
x1996$en.L <- strptime(x1996$ed.L, format ="%Y-%m-%d")$yday + 1

# adding column for phenCt
nf <- aggregate(sd.E ~ cgPlaId, x1996, FUN = length)
names(nf)[names(nf)=="sd.E"] <- "phenCt"
x1996 <- merge(x1996, nf, by.x="cgPlaId", by.y="cgPlaId")

# pedigree data 
ped96 <- ped[ped$cgplaid < 647, ]
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

### get rid of empty levels 
x1996 <- x1996[x1996$site != "", ]
x1996 <- x1996[x1996$site != "Unknown", ]
x1996 <- x1996[x1996$site != "Tower", ] # Only 1 plant flowered
x1996$site <- x1996$site[ , drop=TRUE]
levels(x1996$site)

# adding in row and position data 
x1996 <- merge(x1996, rowpos, by.x = "cgPlaId", by.y = "cgPlaId")

x2013 <- x1996[x1996$year=="2013", ]
x2012 <- x1996[x1996$year=="2012", ]
x2011 <- x1996[x1996$year=="2011", ]
x2010 <- x1996[x1996$year=="2010", ]
x2009 <- x1996[x1996$year=="2009", ]
x2008 <- x1996[x1996$year=="2008", ]
x2007 <- x1996[x1996$year=="2007", ]
x2006 <- x1996[x1996$year=="2006", ]
x2005 <- x1996[x1996$year=="2005", ]

# just getting phenCt
library(plyr)
phenCts <- join(ped96, nf, type = "left")
phenCts <- phenCts[phenCts$site != "", ]
phenCts <- phenCts[phenCts$site != "Unknown", ]
phenCts <- phenCts[phenCts$site != "tower", ] # Only 1 plant flowered
newDF <- subset(phenCts, is.na(phenCts$phenCt))
badPlas <- newDF[,"cgPlaId"]
phenCts <- phenCts[!(phenCts$cgPlaId %in% badPlas), ]
phenCts$site <- phenCts$site[ , drop=TRUE]
levels(phenCts$site)

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

# rank flowering time #####

### first day rank ###
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

x1996 <- do.call(rbind, list(x2005, x2006, x2007, x2008, x2009, x2010, x2011, x2012, x2013))

### DAM rank ####

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

x1996 <- do.call(rbind, list(x2005, x2006, x2007, x2008, x2009, x2010, x2011, x2012, x2013))
str(x1996)
plot(medRank ~ site, x1996)

# Number of flowering plants from each site in each year
aggregate(sn.E ~ year + site, x1996, length)

plot(x2005$rank ~ x2005$site)
plot(x2006$rank ~ x2006$site)
plot(x2007$rank ~ x2007$site)
plot(x2008$rank ~ x2008$site)
plot(x2009$rank ~ x2009$site)
plot(x2010$rank ~ x2010$site)
plot(x2011$rank ~ x2011$site)
plot(x2012$rank ~ x2012$site)
plot(x2013$rank ~ x2013$site)

plot(x2005$medRank ~ x2005$site)
plot(x2006$medRank ~ x2006$site)
plot(x2007$medRank ~ x2007$site)
plot(x2008$medRank ~ x2008$site)
plot(x2009$medRank ~ x2009$site)
plot(x2010$medRank ~ x2010$site)
plot(x2011$medRank ~ x2011$site)
plot(x2012$medRank ~ x2012$site)
plot(x2013$medRank ~ x2013$site)

# Look mean and range of DAM values per plant ####
# remove tower site - only 1 plant 
# remove plants that only flowered once

mm <- x1996[x1996$phenCt != 1, ]
aggregate(medRank ~ cgPlaId, mm, length)
means <- aggregate(medRank ~ cgPlaId, mm, mean)
names(means)[names(means)=="medRank"] <- "meanMedRank"
mm <- merge(mm, means, by.x = "cgPlaId", by.y = "cgPlaId")
plot(meanMedRank ~ site, mm)
aggregate(meanMedRank ~ site, mm, mean)
str(mm)

# adding column for phenCt levels
mm$phenLevel <- "high"
mm$phenLevel[mm$phenCt == 4] <- "mid"
mm$phenLevel[mm$phenCt < 4] <- "low"
mm$phenLevel <- as.factor(mm$phenLevel)
levels(mm$phenLevel)


# visualizing distributions
hist(mm$meanMedRank)
hist(mm[mm$site == "aa", "meanMedRank"], 20)
hist(mm[mm$site == "eri", "meanMedRank"], 20)
hist(mm[mm$site == "lf", "meanMedRank"], 20)
hist(mm[mm$site == "ness", "meanMedRank"], 20)
hist(mm[mm$site == "nwlf", "meanMedRank"], 20)
hist(mm[mm$site == "spp", "meanMedRank"], 20)
hist(mm[mm$site == "sap", "meanMedRank"], 20)
str(mm)

m0 <- lm(meanMedRank ~ 1, mm)
m1 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row + pos + row:pos, mm) # no
m2 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row + pos, mm) # yes
m3 <- lm(meanMedRank ~ site + phenCt + row + pos, mm) # no
m4 <- lm(meanMedRank ~ site + phenCt + site:phenCt + row, mm) # no
m5 <- lm(meanMedRank ~ site + phenCt + site:phenCt + pos, mm) # yes
m6 <- lm(meanMedRank ~ site + phenCt + site:phenCt, mm)
m7 <- lm(meanMedRank ~ site + phenCt + pos, mm)

m8 <- lm(meanMedRank ~ row + pos + row:pos, mm)
m9 <- lm(meanMedRank ~ row + pos, mm)
m10 <- lm(meanMedRank ~ pos, mm)
anova(m9, m10)

anova(m1, m2) # Not significantly different, go with m2
anova(m2, m3) # Yes sig. diff. Interactive site/phenCt term matters. Start with m2
anova(m4, m2) # Yes sig. diff. Pos matters
anova(m5, m2) # Not sig. diff. Row doesn't matter. 
anova(m6, m5) # Pos matters
anova(m7, m5) # interactive term matters

# here's our model  - meanMedRank ~ site + phenCt + site:phenCt + pos

par(mfcol = c(2, 2))
plot(m5)  # looks OK?
par(mfcol = c(1, 1))

plot(meanMedRank ~ phenCt, mm)
a2 <- mm[mm$phenCt==2, ]
a3 <- mm[mm$phenCt==3, ]
a4 <- mm[mm$phenCt==4, ]
a5 <- mm[mm$phenCt==5, ]
a6 <- mm[mm$phenCt==6, ]
a7 <- mm[mm$phenCt==7, ]
a8 <- mm[mm$phenCt==8, ]

plot(meanMedRank ~ site, a2)
plot(meanMedRank ~ site, a3)
plot(meanMedRank ~ site, a4)
plot(meanMedRank ~ site, a5)
plot(meanMedRank ~ site, a6)
plot(meanMedRank ~ site, a7)
plot(meanMedRank ~ site, a8)

coef(m5) # THIS IS CONFUSING?!
summary(m5)

# predicted values
str(mm)

aa <- mm[mm$site=="aa", ]
eri <- mm[mm$site=="eri", ]
lf <- mm[mm$site=="lf", ]
nwlf <- mm[mm$site=="nwlf", ]
ness <- mm[mm$site=="ness", ]
sap <- mm[mm$site=="sap", ]
spp <- mm[mm$site=="spp", ]

range(mm$phenCt)
pv <- 2:8
pv <- rep(pv, 7)

pv2 <- c(rep("ness", 7), rep("lf", 7), rep("sap", 7), 
         rep("spp", 7), rep("aa", 7), rep("nwlf", 7), 
         rep("eri", 7))
pv2
newd <- data.frame(phenCt= pv, site= pv2)
newd
str(newd)
str(mm)
newd$fit <- predict(m6, newd, type = "response") # Using m6 bc doesn't have position
str(newd)

plot(mm$phenCt, mm$meanMedRank)
points(aa$phenCt, aa$meanMedRank, pch = 19, col = "blue")
points(eri$phenCt, eri$meanMedRank, pch = 19, col = "green")
points(ness$phenCt, ness$meanMedRank, pch = 19, col = "purple")
points(lf$phenCt, lf$meanMedRank, pch = 19, col = "red")
points(nwlf$phenCt, nwlf$meanMedRank, pch = 19, col = "gold")
points(spp$phenCt, spp$meanMedRank, pch = 19, col = "pink")
points(sap$phenCt, sap$meanMedRank, pch = 19, col = "gray")

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

# now trying as with phenCt as level 
# 
# m0 <- lm(meanMedRank ~ 1, mm)
# mA <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row + pos + row:pos, mm) # no
# mB <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row + pos, mm) # yes
# mC <- lm(meanMedRank ~ site + phenLevel + row + pos, mm) # no
# mD <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + row, mm) # no
# mE <- lm(meanMedRank ~ site + phenLevel + site:phenLevel + pos, mm) # yes
# mF <- lm(meanMedRank ~ site + phenLevel + site:phenLevel, mm)
# mG <- lm(meanMedRank ~ site + phenLevel + pos, mm)
# 
# anova(mA, mB) # Not significantly different, go with m2
# anova(mB, mC) # Yes sig. diff. Interactive site/phenCt term matters. Start with m2
# anova(mD, mB) # Yes sig. diff. Pos matters
# anova(mE, mB) # Not sig. diff. Row doesn't matter. 
# anova(mF, mE) # Pos matters
# anova(mG, mE) # interactive term matters
# anova(m0, mE) # End up with same model


# testing ranges per site ####

rx <- aggregate(medRank ~ cgPlaId, mm, max)
ry <- aggregate(medRank ~ cgPlaId, mm, min)
rz <- merge(rx, ry, by.x = "cgPlaId", by.y = "cgPlaId")
rz$range <- rz$medRank.x - rz$medRank.y # range of values
mm <- merge(mm, rz, by.x = "cgPlaId", by.y = "cgPlaId")
hist(mm$range, 20, xlab = "Per plant DAM range", main = "")

plot(range ~ site, mm)
aggregate(range ~ site, mm, mean)

nm <- lm(range ~ 1, mm)
am1 <- lm(range ~ site + phenCt + site:phenCt + row + pos + row:pos, mm) # no
am2 <- lm(range ~ site + phenCt + site:phenCt + row + pos, mm) # yes
am3 <- lm(range ~ site + phenCt + row + pos, mm) # 
am4 <- lm(range ~ site + phenCt + site:phenCt + row, mm) # no
am5 <- lm(range ~ site + phenCt + site:phenCt + pos, mm) # yes
am6 <- lm(range ~ site + phenCt + site:phenCt, mm)
am7 <- lm(range ~ site + phenLevel + pos, mm)
am8 <- lm(range ~ site, mm)
am9 <- lm(range ~ pos, mm)

anova(am1, am2)
anova(am2, am3) # site, phenCt matters
anova(am2, am4) # row matters
anova(am2, am5) # pos matters
anova(am2, nm)
anova(am9, nm)

# model = range ~ site + phenCt + site:phenCt + row + pos

aa <- mm[mm$site=="aa", ]
eri <- mm[mm$site=="eri", ]
lf <- mm[mm$site=="lf", ]
nwlf <- mm[mm$site=="nwlf", ]
ness <- mm[mm$site=="ness", ]
sap <- mm[mm$site=="sap", ]
spp <- mm[mm$site=="spp", ]

range(mm$phenCt)
pv <- 2:8
pv <- rep(pv, 7)

pv2 <- c(rep("ness", 7), rep("lf", 7), rep("sap", 7), 
         rep("spp", 7), rep("aa", 7), rep("nwlf", 7), 
         rep("eri", 7))
pv2
newd <- data.frame(phenCt= pv, site= pv2)
newd
str(newd)
str(mm)
newd$fit <- predict(am6, newd, type = "response")
str(newd)

plot(mm$phenCt, mm$range)
points(aa$phenCt, aa$range, pch = 19, col = "blue")
points(eri$phenCt, eri$range, pch = 19, col = "green")
points(ness$phenCt, ness$range, pch = 19, col = "purple")
points(lf$phenCt, lf$range, pch = 19, col = "red")
points(nwlf$phenCt, nwlf$range, pch = 19, col = "gold")
points(spp$phenCt, spp$range, pch = 19, col = "pink")
points(sap$phenCt, sap$range, pch = 19, col = "gray")

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


par(mfcol = c(2,2))
plot(m2)
par(mfcol = c(1,1))

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

aggregate(cgPlaId ~ site, mm, length)




