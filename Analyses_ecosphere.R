rm(list=ls()); gc()

### Packages
require(vegan)
require(cooccur)
require(grid)
require(ggpubr)
require(netassoc)
require(factoextra)
require(labdsv)
require(rio)
require(data.table)
require(tidyverse)
require(pheatmap)
require(ade4)
require(adegraphics)
require(multipatt)
require(gridExtra)

### Set working directory
path <- "Z:/intranet/espaces_individuels_lebco/mchevalier/Underprocessed papers/Bunyeth_co-occurrence/Analyses/"
setwd(path)

### Get data
data <- rio::import("data_for_analyses.txt")
data$period <- as.factor(data$period)
levels(data$period) <- c("High", "Low", "Receding", "Rising")

### Load information on migratory status
sp.features <- import("species_bchan.xlsx")
sp.features$Species <- gsub(" ", ".", sp.features$Species)
names.estuarine <- sp.features[which(sp.features$`Migratory pattern`=="Estuarine"), "Species"]
names.resident <- sp.features[which(sp.features$`Migratory pattern`=="floodplain resident"), "Species"]
names.lateral <- sp.features[which(sp.features$`Migratory pattern`=="Lateral"), "Species"]
names.longitudinal <- sp.features[which(sp.features$`Migratory pattern`=="longitudinal"), "Species"]

### Remove estuarine species
data <- data[,-which(colnames(data) %in% names.estuarine)]

### Keep species matching with the two databases
data <- data[,c(1:5, which(colnames(data) %in% c(names.resident, names.lateral, names.longitudinal)))]
Nsp <- ncol(data[,-c(1:5)])

#######################
### Analysis on abundances
#######################

data.ab <- data
data.tot.ab <- apply(data.ab[,-c(1:5)], 1, sum)
data.tot.ab.resident <- apply(data.ab[,which(colnames(data.ab) %in% names.resident)], 1, sum)
data.tot.ab.lateral <- apply(data.ab[,which(colnames(data.ab) %in% names.lateral)], 1, sum)
data.tot.ab.longitudinal <- apply(data.ab[,which(colnames(data.ab) %in% names.longitudinal)], 1, sum)
df.ab <- as.data.frame(cbind(data.tot.ab, data.tot.ab.resident, data.tot.ab.lateral, data.tot.ab.longitudinal))
colnames(df.ab) <- c("Total", "Resident", "Lateral", "Longitudinal")
df.ab <- cbind(data.ab[,2:5], df.ab)

#----------------
### Compare abundances between periods for each species
#----------------

# test differences for each species between periods
pval <- numeric()
for(i in 1:Nsp){
  ks <- kruskal.test(data.ab[,(i+5)], data.ab$period)
  pval[i] <- ks$p.value
}
padj <- p.adjust(pval)
names(padj) <- colnames(data.ab[,-c(1:5)])
posi <- which(padj<0.05)

# Visualise periodic changes in abundances
df.ab.melt <- melt(data.ab[,-c(1,3,5)], id.vars=c("site", "period"))
agg <- aggregate(df.ab.melt$value, by=list(df.ab.melt$period, df.ab.melt$variable), FUN=mean)  
mat <- match(agg$Group.2, names(padj[posi]))
agg$Group.3 <- 0
agg$Group.3[which(!is.na(mat))] <- 1
colnames(agg)[4] <- "Test"
agg$Test <- factor(agg$Test)
plotS1 <- ggplot(agg, aes(x=Group.1, y=Group.2, size=log1p(x), col=Test)) +
  geom_point(alpha=.2) +
  theme_bw()+
  scale_size_area(breaks = seq(0,7,1))+
  labs(y="", x="")

tiff("FigS1_Bubble_abundance_signif_change_period.tiff", res=250, width=1800, height=2000)
print(plotS1)
dev.off()

### Plot total abundance across periods for each migratory status
df.ab.melt <- melt(df.ab[,-5], id.vars=c("site", "date", "period", "season"))
plotS2 <- ggplot(df.ab.melt, aes(x=variable, y=log(value+1), col=period))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5))+
  theme_bw()+
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a"))+
  labs(x="", y="Log(abundance+1)")
tiff("FigS2_Abundance_migratory_period.tiff", res=250, width=1300, height=1000)
print(plotS2)
dev.off()

### Multivariate analysis (FIGURE 3)
nmds.ab <- metaMDS(data.ab[,-c(1:5)], distance = "bray")
Permanova.ab <- adonis2(data.ab[,-c(1:5)] ~ data.ab$period * data.ab$site, permutations=999)
coord.ab <- nmds.ab$points[,c(1,2)]
species.scores <- as.data.frame(vegan::scores(nmds.ab, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$Species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores <- merge(species.scores, sp.features[,c(2,7)])
colnames(species.scores)[4] <- "Mig"

#-- Period
g1 <- s.class(coord.ab, fac=factor(data.ab$period), col=c("#1b9e77","#d95f02","#7570b3","#e7298a"), starSize = 0, ppoints.cex=0)

#-- sites
g2 <- s.class(coord.ab, fac=factor(data.ab$site), col=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"), starSize = 0, ppoints.cex=0)

#-- Migratory status
g3 <- s.class(species.scores[,2:3], fac=factor(species.scores$Mig),col=c("#bebada","#fb8072","#b3de69"), 
        ppoints.cex=0, starSize = 0, samelimits = FALSE)

#-- Combined figure
tiff("nmds_abundance_per_site_mig.tiff", res=250, width=2000, height=800)
ADEgS(adeglist = list(g1, g2, g3))
dev.off()

### Site-wise analysis
tiff("nmds_abundance_community_site-wise.tiff", res=250, width=2000, height=800)
g4 <- s.class(coord.ab, fac=factor(data.ab$period), col=c("#1b9e77","#d95f02","#7570b3","#e7298a"), 
              ppoints.cex=0, starSize = 0, facets=data.ab$site, samelimits = FALSE, pfacets.cex=2)
dev.off()

#######################
### Analysis on occurences
#######################

data.pa <- data.ab
data.pa[,-c(1:5)] <- ifelse(data.pa[,-c(1:5)]>0, 1, 0)
data.tot.pa <- apply(data.pa[,-c(1:5)], 1, sum)
data.tot.pa.resident <- apply(data.pa[,which(colnames(data.pa) %in% names.resident)], 1, sum)
data.tot.pa.lateral <- apply(data.pa[,which(colnames(data.pa) %in% names.lateral)], 1, sum)
data.tot.pa.longitudinal <- apply(data.pa[,which(colnames(data.pa) %in% names.longitudinal)], 1, sum)
df.pa <- as.data.frame(cbind(data.tot.pa, data.tot.pa.resident, data.tot.pa.lateral, data.tot.pa.longitudinal))
colnames(df.pa) <- c("Total", "Resident", "Lateral", "Longitudinal")
df.pa <- cbind(data.pa[,2:5], df.pa)
# df.pa$date.id <- as.numeric(df.pa$date)

# ### Plot richness across periods and Sites
# ggplot(df.pa, aes(x=period, y=Total, col=site))+
#   stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
#   stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5))+
#   theme_bw()+
#   labs(x="", y="Species richness")

### Plot richness across periods for each migratory status (FIGURE 2A)
df.pa.melt <- melt(df.pa[,-5], id.vars=c("site", "date", "period", "season"))
plot1 <- ggplot(df.pa.melt, aes(x=variable, y=value, fill=period))+
  geom_boxplot(position=position_dodge(width=.5))+
  # stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) +
  stat_summary(fun.y = mean, geom ="point", size = 2, show.legend = FALSE, position=position_dodge(width=.5), col="grey80")+
  theme_bw()+
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a"))+
  labs(x="", y="Species richness")

tiff("Final_figures/Fig2.tiff", res=300, width=2200, height=1800)
print(plot1)
dev.off()

### Plot relative frequency of occurrence for a set of longitudinal species (FIGURE 2B)
tmp <- indval(data.pa[,-c(1:5)], clustering=data.pa$period)
tmp <- as.data.frame(tmp$relfrq)
tmp$species <- rownames(tmp)

tmp_long <- tmp %>%
  pivot_longer(cols = -species,
               names_to = "variable",
               values_to = "value")

mat <- xtabs(value ~ variable + species, data = tmp_long)
Heat0 <- pheatmap(mat, angle_col = 90, border_color = NA,  fontsize_col = 6)
Heat0_grob <- Heat0$gtable 

tiff("Final_figures/Fig3.tiff", res=350, width=2600, height=2000)
pheatmap(mat, angle_col = 90, border_color = NA,  fontsize_col = 6)
dev.off()

# tiff("Fig2_richness_heatmap.tiff", res=250, width=3800, height=1800)
# # grid.arrange(plot1, plot2, ncol=2)
# ggarrange(plot1, Heat0_grob, ncol=2, nrow=1, labels = c("(A)", "(B)"))
# dev.off()

# ### bubbble plot showing the relative frequency of occurrence for all species
# tmp <- indval(data.pa[,-c(1:5)], clustering=data.pa$period)
# tmp <- tmp$relfrq
# tmp$species <- rownames(tmp)
# tmp <- melt(tmp, id.vars="species")
# ggplot(tmp, aes(x=variable, y=species, size=sqrt(value))) +
#   geom_point(alpha=.2) +
#   theme_bw()

# ### Plot temporal changes in abundance of white, grey and black fishes
# ggplot(df.pa.melt, aes(x=date, y=value))+ #, col=period
#   geom_point(alpha=.2, size=.4)+
#   geom_smooth(method='gam')+
#   theme_bw()+
#   theme(legend.position = "none")+
#   facet_wrap(~variable, scale="free_y")+
#   labs(y="Species richness", x="Time")

### Plot proportion of each migratory species across sites
tmp <- df.pa
tmp[,6:8] <- apply(tmp[,6:8], 2, function(x)(x/tmp$Total)*100)
tmp.melt <- melt(tmp[,-5], id.vars=c("site", "date", "period", "season"))
plotS3 <- ggplot(tmp.melt, aes(x=variable, y=value, col=site))+
  # geom_boxplot(position=position_dodge(width=.5))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) +
  stat_summary(fun = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5))+
  theme_bw()+
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")) +
  labs(x="", y="Proportion of species")

tiff("FigS3_prop_spMig_site.tiff", res=250, width=1300, height=1000)
print(plotS3)
dev.off()

# ### Bubble plot (Which species are specific to which season?)
# df.pa.melt <- melt(data.pa[,-c(1,3,5)], id.vars=c("site", "period"))
# agg <- aggregate(df.pa.melt$value, by=list(df.pa.melt$period, df.pa.melt$variable), FUN=mean)  
# ggplot(agg, aes(x=Group.1, y=Group.2, size=x)) +
#   geom_point(alpha=.2) +
#   theme_bw()+
#   scale_size_area(breaks = seq(0,7,1))

# ### Multivariate analysis
# nmds.pa <- metaMDS(data.pa[,-c(1:5)], distance = "jaccard")
# Permanova.pa <- adonis2(data.pa[,-c(1:5)] ~ data.pa$period + data.pa$site + data.pa$season, permutations=500)
# coord.pa <- nmds.pa$points[,c(1,2)]
# s.class(coord.pa, fac=factor(data.pa$site), col=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"), possub="topleft", sub="(a)", ppoints.cex=0.2)
# s.class(coord.pa, fac=factor(data.pa$period), col=c("#1b9e77","#d95f02","#7570b3","#e7298a"), possub="topleft", sub="(b)", ppoints.cex=0, starSize = 0)
# s.class(coord.pa, fac=factor(data.pa$period), col=c("#1b9e77","#d95f02","#7570b3","#e7298a"), ppoints.cex=0, starSize = 0, facets=data.pa$site, samelimits = FALSE)

#######################
### spatial synchrony on population abundances - variation across periods
#######################

synch <- NULL
df.ab <- data.ab[,-c(1:5)]
periods <- levels(factor(data.ab$period))
seasons <- levels(factor(data.ab$season))
df.ab <- data.ab[,c(1:5, which(colnames(data.ab) %in% sp.features$Species))]
Nsp <- ncol(df.ab) - 5
sites <- levels(factor(df.ab$site))

for(i in 1:length(periods)){
  
  df.period <- df.ab[which(df.ab$period == periods[i]),]
  
  # find dates in common between all TS
  site.1 <- as.character(df.period[which(df.period$site == "bb"),"date"])
  site.2 <- as.character(df.period[which(df.period$site == "kc"),"date"])
  mat <- intersect(site.1, site.2)
  for(z in 3:length(sites)){
    site.z <- as.character(df.period[which(df.period$site == sites[z]),"date"])
    mat <- intersect(mat, site.z)
  }
  
  # Extract site id for the common dates
  df.period.common <- df.period[which(df.period$date %in% mat),]
  id.bb <- which(df.period.common$site == "bb")
  id.kc <- which(df.period.common$site == "kc")
  id.kt <- which(df.period.common$site == "kt")
  id.ps <- which(df.period.common$site == "ps")
  id.sr <- which(df.period.common$site == "sr")
  
  # Compute synchrony
  for(j in 1:Nsp){
    
    focal.sp <- cbind(df.period.common[id.bb,(j+5)], 
                      df.period.common[id.kc,(j+5)],
                      df.period.common[id.kt,(j+5)],
                      df.period.common[id.ps,(j+5)],
                      df.period.common[id.sr,(j+5)])
    colnames(focal.sp) <- c("bb", "kc", "kt", "ps", "sr")
    
    id.pres <- which(apply(focal.sp, 2, function(x)length(which(x!=0))/length(x)) > 0.5) # Select TS with 50% completeness
    mig.status <- sp.features[which(sp.features$Species == colnames(df.period.common)[(j+5)]), "Migratory pattern"]
    
    if(length(id.pres) > 1){
      focal.sp <- focal.sp[,id.pres]
      correl <- cor(focal.sp, method="spearman")
      
      ### If we want to averaghe correlations at the spaecies scale
      # avg.correl <- mean(correl[upper.tri(correl, diag = FALSE)])
      # out <- c(periods[i], colnames(data.ab)[(j+5)], length(id.pres), mig.status, avg.correl)
      
      ### If we want to keep information on synchrony between populations
      correl <- correl[upper.tri(correl, diag = FALSE)]
      out <- cbind(periods[i], colnames(data.ab)[(j+5)], length(id.pres), mig.status, correl)
    } else {
      out <- c(periods[i], colnames(data.ab)[(j+5)], length(id.pres), mig.status, NA)
    }
    
    synch <- rbind(synch, out)
  }
  
}

synch <- as.data.frame(synch)
colnames(synch) <- c("Periods", "Species", "N.sites", "mig.status", "Synchrony")
synch$Synchrony <- as.numeric(synch$Synchrony)
synch$N.sites <- as.numeric(synch$N.sites)
synch$Periods <- factor(synch$Periods)
synch$mig.status <- factor(synch$mig.status)
synch <- na.omit(synch)
table(synch$Periods, synch$mig.status)
synch <- droplevels(synch)
table(synch$mig.status)
tapply(synch$Synchrony, synch$Species, mean)

### Draw heatmaps
require(pheatmap)
synch$Periods <- factor(synch$Periods, levels = c("Low", "Rising", "High", "Receding"))
mat <- aggregate(Synchrony ~ Species + Periods, data = synch, FUN = mean)
mat <- xtabs(Synchrony ~ Species + Periods, data = mat)
lim <- max(abs(mat), na.rm = TRUE)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
my_breaks <- seq(-lim, lim, length.out = 101)

tiff("Heatmap_synch_periods.tiff", res=250, width=2000, height=1800)
Heat1 <- pheatmap(mat, color = my_colors, breaks = my_breaks, cluster_rows = FALSE, border_color = NA,
                  cluster_cols = FALSE, na_col = "grey80",  angle_col = 0)
dev.off()
Heat1_grob <- Heat1$gtable 

# kruskal.test(synch$Synchrony ~ synch$mig.status)
anova(lm(synch$Synchrony ~ synch$Periods * synch$mig.status, weight=synch$N.sites))
# require(car)
# leveneTest(synch$Synchrony ~ synch$Periods * synch$mig.status)

plot1 <- ggplot(synch, aes(x=mig.status, y=Synchrony, col=Periods))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
  theme_bw()+
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a"))+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(x="", y="Within species synchrony", col="Migratory guild")

#######################
### Site-wise species synchrony on population abundances 
#######################

synch.sp <- NULL
periods <- levels(factor(data.ab$period))
sites <- levels(factor(data.ab$site))
df.ab <- data.ab[,c(1:5, which(colnames(data.ab) %in% sp.features$Species))]

for(i in 1:length(periods)){
  
  df.period <- df.ab[which(df.ab$period == periods[i]),]
  
  for(j in 1:length(sites)){
    
    df.site <- df.period[which(df.period$site == sites[j]),]
    # df.site <- df.period
    id.pres <- which(apply(df.site[,-c(1:5)], 2, function(x)length(which(x!=0))/length(x)) > 0.5) # Select TS with 50% completeness
    sp.list <- names(id.pres)
    
    for(k1 in 1:(length(sp.list)-1)){
      
      for(k2 in (k1+1):length(sp.list)){
        
        mig.status.sp1 <- mig.status.sp2 <- cor.sp <- id.mig <- NA
        
        focal.sp1 <- df.site[,sp.list[k1]]
        focal.sp2 <- df.site[,sp.list[k2]]
        cor.sp <- cor(focal.sp1, focal.sp2, method="spearman")
        mig.status.sp1 <- sp.features[which(sp.features$Species == sp.list[k1]), "Migratory pattern"]
        mig.status.sp2 <- sp.features[which(sp.features$Species == sp.list[k2]), "Migratory pattern"]
        id.mig <- mig.status.sp1 == mig.status.sp2 # define if they have the same migratory status or not (TRUE= same, FALSE=different)
        
        out <- c(periods[i], sites[j], sp.list[k1], sp.list[k2], cor.sp, mig.status.sp1, mig.status.sp2, id.mig)
        synch.sp <- rbind(synch.sp, out)
      }
      
    }
    
  }
  
}

TMP <- sp.features[,c(2,3,4,5,6,7)]
TMP <- TMP[which(TMP$Species %in% colnames(data.ab)),]
write.table(TMP[,-1], file="Table_S1.txt")

synch.sp <- as.data.frame(synch.sp)
colnames(synch.sp) <- c("Periods", "Sites", "Sp1", "Sp2", "Synch", "mig.status.sp1", "mig.status.sp2", "ID.mig")
synch.sp$Synch <- as.numeric(synch.sp$Synch)
synch.sp$Periods <- factor(synch.sp$Periods)
synch.sp$Sites <- factor(synch.sp$Sites)
synch.sp$ID.mig <- factor(synch.sp$ID.mig)
levels(synch.sp$ID.mig) <- c("no","yes")

range(synch.sp$Synch)
mean(synch.sp$Synch)
sd(synch.sp$Synch)

plotS4 <- ggplot(synch.sp, aes(x=Sites, col=Periods, y=Synch))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
  theme_bw()+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(x="", y="Among species synchrony", col="Periods")

tiff("FigS4_BtwSynch_site_mig.tiff", res=250, width=1300, height=1000)
print(plotS4)
dev.off()

# table(synch.sp$ID.mig, synch.sp$Periods)
anova(lm(synch.sp$Synch ~ synch.sp$Periods * synch.sp$ID.mig))

lim <- max(abs(synch.sp$Synch), na.rm = TRUE)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
my_breaks <- seq(-lim, lim, length.out = 101)
per.lev <- c("Low", "Rising", "High", "Receding")
plots <- list()
for(i in 1:4){
  tmp <- synch.sp[which(synch.sp$Periods == per.lev[i]),]
  mat <- aggregate(Synch ~ Sp1 + Sp2, data = tmp, FUN = mean)
  mat <- xtabs(Synch ~ Sp1 + Sp2, data = mat)
  plots[[i]] <- pheatmap(mat, color = my_colors, breaks = my_breaks, cluster_rows = FALSE, 
                         cluster_cols = FALSE, na_col = "grey80", main=per.lev[i])
}

tiff("Heatmap_btw_synch_periods.tiff", res=250, width=4000, height=3500)
Heat2 <- grid.arrange(grobs = lapply(plots, function(x) x$gtable), ncol = 2)   # choose number of columns
dev.off()


plot2 <- ggplot(synch.sp, aes(x=Periods, col=ID.mig, y=Synch))+
  # geom_boxplot(position=position_dodge(width=.5))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
  theme_bw()+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_manual(values = c("#dfc27d","#80cdc1"))+
  labs(x="", y="Among species synchrony", col="Same guild")#+
  # facet_wrap(~Sites)

plot2bis <- ggplot(synch.sp[which(synch.sp$ID.mig=="yes"),], aes(x=Periods, col=mig.status.sp1, y=Synch))+
  # geom_boxplot(position=position_dodge(width=.5))+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
  theme_bw()+
  labs(x="", y="Among species synchrony", col="Migratory guild")+
  scale_color_manual(values = c("#bebada","#fb8072","#b3de69"))+
  geom_hline(yintercept = 0, linetype="dashed")
# facet_wrap(~Sites)

tiff("Synchrony.tiff", res=250, width=3000, height=2000)
# grid.arrange(plot1, plot1bis, plot2, plot2bis, ncol=2)
ggarrange(plot1, Heat1_grob, plot2, plot2bis, ncol=2, nrow=2, labels = c("(A)", "(B)", "(C)", "(D)"))
dev.off()

##################################
### Cooccurence analysis
##################################

#--------- Whole data
z <- cooccur(t(data.pa[,-c(1:5)]), type = "spp_site", thresh = T, spp_names=T)
prob.table(z)
profile <- pair.profile(z)
exp.exp <- obs.v.exp(z)
plot(z)
prop.pos <- (z$positive/z$pairs)*100 # proportion of positive associations
prop.neg <- (z$negative/z$pairs)*100 # proportion of negative associations
prop.rand <- 100-(prop.pos+prop.neg) # proportion of random associations

### Do associations vary depending on whether you are in the same guild or not?
tmp.sp1 <- match(z$results$sp1_name, sp.features$Species)
tmp.sp2 <- match(z$results$sp2_name, sp.features$Species)
mig.status.sp1 <- sp.features[tmp.sp1, "Migratory pattern"]
mig.status.sp2 <- sp.features[tmp.sp2, "Migratory pattern"]
id.mig <- mig.status.sp1 == mig.status.sp2
A <- z$results
A$Mig <- id.mig
A$assoc <- 0 # random
A[which(A$p_lt < 0.05), "assoc"] <- (-1) # negative
A[which(A$p_gt < 0.05), "assoc"] <- 1 # positive
tmp <- table(A$assoc, A$Mig)
apply(tmp,2,function(x) x/sum(x)) # NOPE it does not vary

# ### Possibility to use heatmaps to plot results
z <- cooccur(t(data.pa[,-c(1:5)]), type = "spp_site", thresh = F, spp_names=T, only_effects=TRUE, eff_standard=TRUE, eff_matrix=TRUE)
heatmap(as.matrix(z), symm = TRUE, cexRow=.5, cexCol=.5)

# ### Possibility to use hierarchical clustering to visualize clusters
# res <- hcut(z, k = 3, stand = TRUE)
# fviz_dend(res, rect = TRUE, cex = 0.8, type="rectangle", k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))  

##########################
#--------- For each period (for a similar analyses across sites use the script "co-occur_site")
##########################

lev.period <- levels(factor(data.pa$period))
prop.pos <- prop.neg <- prop.rand <- c()
list.cooc <- list()
list.profile <- list()
sp.names <- colnames(data.pa[,-c(1:5)])

for(i in 1:length(lev.period)){
  
  ### Co-occurrence analysis
  z <- cooccur(t(data.pa[which(data.pa$period == lev.period[i]),-c(1:5)]), type = "spp_site", thresh = F, spp_names = T, prob="comb")
  prop.pos <- c(prop.pos, (z$positive/z$pairs)*100)
  prop.neg <- c(prop.neg, (z$negative/z$pairs)*100)
  prop.rand <- c(prop.rand, 100-(prop.pos[i]+prop.neg[i]))
  list.cooc[[i]] <- z
  list.profile[[i]] <- pair.profile(z)
  
  ### Matrix looking at changes in species associations between periods
  cooc <- list.cooc[[i]]$results
  comat_pos <- comat_neg <- matrix(0, nrow = list.cooc[[i]]$species, ncol = list.cooc[[i]]$species)
  colnames(comat_pos) <- colnames(comat_neg) <- sp.names
  rownames(comat_pos) <- rownames(comat_neg) <- sp.names
  diag(comat_pos) <- diag(comat_neg) <- NA
  for (j in 1:nrow(cooc)) {
    comat_pos[which(rownames(comat_pos)==cooc[j, "sp1_name"]), which(colnames(comat_pos)==cooc[j, "sp2_name"])] <- cooc[j, "p_gt"]
    comat_pos[which(rownames(comat_pos)==cooc[j, "sp2_name"]), which(colnames(comat_pos)==cooc[j, "sp1_name"])] <- cooc[j, "p_gt"]
    comat_neg[which(rownames(comat_neg)==cooc[j, "sp1_name"]), which(colnames(comat_neg)==cooc[j, "sp2_name"])] <- cooc[j, "p_lt"]
    comat_neg[which(rownames(comat_neg)==cooc[j, "sp2_name"]), which(colnames(comat_neg)==cooc[j, "sp1_name"])] <- cooc[j, "p_lt"]
  }
  comat <- ifelse(comat_pos <= 0.05, 1, 0) + ifelse(comat_neg <= 0.05, -1, 0)
  assign(paste0("comat.", lev.period[i]), comat)
  
}
names(list.cooc) <- names(list.profile) <- lev.period

#============= Plot changes in species associations

### Positive associations
df.assoc.pos <- as.data.frame(cbind(apply(comat.High, 2, function(x)length(which(x==1))/(length(x)-1)), # -1 to remove the NA in the diag 
                                    apply(comat.Low, 2, function(x)length(which(x==1))/(length(x)-1)),
                                    apply(comat.Receding, 2, function(x)length(which(x==1))/(length(x)-1)),
                                    apply(comat.Rising, 2, function(x)length(which(x==1))/(length(x)-1))))
colnames(df.assoc.pos) <- lev.period
df.assoc.pos$Species <- colnames(data.ab[,-c(1:5)])
df.assoc.pos <- df.assoc.pos[which(df.assoc.pos$Species %in% sp.features$Species),]
df.assoc.pos$mig <- sp.features[which(sp.features$Species %in% df.assoc.pos$Species),"Migratory pattern"]

### Random associations
df.assoc.rand <- as.data.frame(cbind(apply(comat.High, 2, function(x)length(which(x==0))/(length(x)-1)),
                                     apply(comat.Low, 2, function(x)length(which(x==0))/(length(x)-1)),
                                     apply(comat.Receding, 2, function(x)length(which(x==0))/(length(x)-1)),
                                     apply(comat.Rising, 2, function(x)length(which(x==0))/(length(x)-1))))
colnames(df.assoc.rand) <- lev.period
df.assoc.rand$Species <- colnames(data.ab[,-c(1:5)])
df.assoc.rand <- df.assoc.rand[which(df.assoc.rand$Species %in% sp.features$Species),]
df.assoc.rand$mig <- sp.features[which(sp.features$Species %in% df.assoc.rand$Species),"Migratory pattern"]

### Negative associations
df.assoc.neg <- as.data.frame(cbind(apply(comat.High, 2, function(x)length(which(x==(-1)))/(length(x)-1)),
                                    apply(comat.Low, 2, function(x)length(which(x==(-1)))/(length(x)-1)),
                                    apply(comat.Receding, 2, function(x)length(which(x==(-1)))/(length(x)-1)),
                                    apply(comat.Rising, 2, function(x)length(which(x==(-1)))/(length(x)-1))))
colnames(df.assoc.neg) <- lev.period
df.assoc.neg$Species <- colnames(data.ab[,-c(1:5)])
df.assoc.neg <- df.assoc.neg[which(df.assoc.neg$Species %in% sp.features$Species),]
df.assoc.neg$mig <- sp.features[which(sp.features$Species %in% df.assoc.neg$Species),"Migratory pattern"]

########################################

# names.resident <- sp.features[which(sp.features$`Migratory pattern`=="floodplain resident"), "Species"]
# names.lateral <- sp.features[which(sp.features$`Migratory pattern`=="Lateral"), "Species"]
# names.longitudinal <- sp.features[which(sp.features$`Migratory pattern`=="longitudinal"), "Species"]
# colnames(comat.High) <- colnames(comat.Low) <- colnames(comat.Receding) <- colnames(comat.Rising) <- colnames(data.ab[,-c(1:5)])
# rownames(comat.High) <- rownames(comat.Low) <- rownames(comat.Receding) <- rownames(comat.Rising) <- colnames(data.ab[,-c(1:5)])
# pos.res <- which(rownames(comat.High) %in% names.resident)
# pos.lat <- which(rownames(comat.High) %in% names.lateral)
# pos.long <- which(rownames(comat.High) %in% names.longitudinal)
# 
# ### Plot positive associations between periods with a 1-1 line
# count <- 0
# for(i in 1:(ncol(df.assoc.pos[,-c(5,6)])-1)){
#   for(j in (i+1):ncol(df.assoc.pos[,-c(5,6)])){
#       tmp.df <- df.assoc.pos[,c(i,j,6)]
#       count <- count + 1
#       tmp.plot <- ggplot(data=tmp.df, aes_string(x=eval(colnames(tmp.df)[1]), y=eval(colnames(tmp.df)[2])))+ 
#         geom_point(data=tmp.df, aes(col=mig))+
#         geom_abline(intercept = 0, slope = 1)+
#         labs(x=colnames(tmp.df)[1], y=colnames(tmp.df)[2])+
#         lims(x=c(0, max(df.assoc.pos[,-c(5,6)])), y=c(0, max(df.assoc.pos[,-c(5,6)])))+
#         theme_bw()
#       assign(paste0("plot", count), tmp.plot)
#   }
# }
# ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

### Plot positive, negative and random associations for each migratory guild and period
df.assoc.neg.melt <- melt(df.assoc.neg, id.vars=c("Species", "mig"))
df.assoc.neg.melt$Assoc <- "Negative"
df.assoc.pos.melt <- melt(df.assoc.pos, id.vars=c("Species", "mig"))
df.assoc.pos.melt$Assoc <- "Positive"
df.assoc.rand.melt <- melt(df.assoc.rand, id.vars=c("Species", "mig"))
df.assoc.rand.melt$Assoc <- "Random"
df.assoc.glob <- rbind(df.assoc.neg.melt, df.assoc.pos.melt, df.assoc.rand.melt) # 
df.assoc.glob$mig <- "Global"
df.tot <- rbind(df.assoc.neg.melt, df.assoc.pos.melt, df.assoc.rand.melt, df.assoc.glob)

df.tot$mig <- factor(df.tot$mig, levels = c("Global", "floodplain resident", "Lateral", "longitudinal"))
df.tot$variable <- factor(df.tot$variable, levels = c("High", "Receding", "Low", "Rising"))

plot.glob <- ggplot(data=df.tot, aes(x=variable, y=value*100, col=Assoc, group=Assoc))+ 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
  stat_summary(fun.y = mean, geom ="line", size = 1,show.legend = FALSE, position=position_dodge(width=.5)) + 
  theme_bw()+
  labs(y="Proportion of associations (%)", x="")+
  facet_wrap(~mig)+
  scale_color_manual(values = c("#33a02c", "#fdbf6f", "#cab2d6"))+
  labs(col="Associations")#+
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tiff("Prop_Assoc.tiff", res=250, width=2000, height=1600)
print(plot.glob)
dev.off()

# ### Plot differences in associations between periods for each type of associations
# df.dif <- NULL
# for(i in 1:(ncol(df.assoc.pos[,-c(5,6)])-1)){
#   for(j in (i+1):ncol(df.assoc.pos[,-c(5,6)])){
#     dif.pos <- df.assoc.pos[,i] - df.assoc.pos[,j]
#     dif.neg <- df.assoc.neg[,i] - df.assoc.neg[,j]
#     dif.rand <- df.assoc.rand[,i] - df.assoc.rand[,j]
#     per <- paste(colnames(df.assoc.pos)[i], colnames(df.assoc.pos)[j], sep="-")
#     df.dif <- rbind(df.dif, cbind(per, df.assoc.pos[,5:6], dif.pos, dif.neg, dif.rand))
#   }
# }
# df.dif <- as.data.frame(df.dif)
# colnames(df.dif)[4:6] <- c("Positive", "Negative", "Random") 
# df.dif <- data.table::melt(df.dif, id.vars=c("per", "Species", "mig"))
# 
# ### Plot for all pairs of season
# tmp.plot <- ggplot(data=df.dif, aes(x=per, y=value, col=mig))+ 
#   stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
#   stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
#   theme_bw()+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   labs(y="Change in associations", x="")+
#   facet_wrap(~variable)+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# tiff("Change_associations.tiff", res=250, width=3000, height=1000)
# print(tmp.plot)
# dev.off()

# ### Plot only between consecutive seasons (high -> receding -> low -> rising)
# df.dif <- df.dif[which(df.dif$per %in% c("High-Receding", "Low-Receding", "Low-Rising", "High-Rising")),]
# df.dif[which(df.dif$per == "Low-Receding"),"value"] <- df.dif[which(df.dif$per == "Low-Receding"),"value"]*(-1)
# df.dif[which(df.dif$per == "High-Rising"),"value"] <- df.dif[which(df.dif$per == "High-Rising"),"value"]*(-1)
# df.dif$per <- factor(df.dif$per)
# levels(df.dif$per) <- c("High-Receding", "Rising-High", "Receding-Low", "Low-Rising")
# df.dif$per <- factor(df.dif$per, levels = c("High-Receding", "Receding-Low", "Low-Rising", "Rising-High"))
# tmp.plot <- ggplot(data=df.dif, aes(x=per, y=value, col=mig))+ 
#   stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
#   stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
#   stat_summary(aes(group=mig), fun.y = mean, geom ="line", show.legend = FALSE, position=position_dodge(width=.5)) + 
#   theme_bw()+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   labs(y="Change in associations", x="")+
#   facet_wrap(~variable)+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# tiff("Change_associations_consecutive_periods.tiff", res=250, width=3000, height=1000)
# print(tmp.plot)
# dev.off()
# 
# ### No discrimination depending on migratory status
# tmp.plot <- ggplot(data=df.dif, aes(x=per, y=value))+ 
#   stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0, size=0.8, position=position_dodge(width=.5)) + 
#   stat_summary(fun.y = mean, geom ="point", size = 2,show.legend = FALSE, position=position_dodge(width=.5)) + 
#   theme_bw()+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   labs(y="Change in associations", x="")+
#   facet_wrap(~variable)+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# tiff("Change_associations_consecutive_periods_global.tiff", res=250, width=3000, height=1000)
# print(tmp.plot)
# dev.off()


#============= Plot changes in species profiles
df.profile <- NULL
name <- names(list.profile)
for(i in 1:4){
  tmp.df <- list.profile[[i]]$data
  tmp.df$period <- name[i]
  df.profile <- rbind(df.profile, tmp.df)
}
df.profile <- as.data.frame(df.profile)
colnames(df.profile)[2] <- "Assoc"
levels(df.profile$Assoc) <- c("Positive", "Negative", "Random")

### Each species
tmp <- df.profile[-which(df.profile$sppname=="All Species"),]
tmp.plot <- ggplot(tmp, aes(x=value, y=sppname, fill=Assoc))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a"))+
  theme_bw()+
  facet_wrap(~period, nrow=1)+
  scale_x_continuous(expand = c(0, 0))+
  labs(x="", y="Percent of pairings", fill="Type of association")

tiff("Final_figures/Fig7.tiff", res=400, width=4800, height=4000)
print(tmp.plot)
dev.off()


### All species
tmp <- df.profile[which(df.profile$sppname=="All Species"),]
tmp.plot <- ggplot(tmp, aes(x=period, y=value, fill=variable))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a"))+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0))+labs(x="")


