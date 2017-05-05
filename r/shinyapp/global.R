# CRAN libraries
library(shiny) # app
library(ggplot2) # plotting
library(plotly) # interactive plots over ggplots
library(reshape2)
library(dplyr)
library(RPostgreSQL)
library(DT)
library(grofit) # calculation of growth slopes
# Bioconductory libraries
library(ggbio) #genome plots

# library(ensembldb) # functions to create and use transcript centric annotation databases/packages
# library(EnsDb.Hsapiens.v75) # ensembl data

# --------------------------------------------------------------------------------
# Should be loaded from DB
# Things annotated as %DB% and /%DB% require updating after NGS data is in the database
# %DB%
# Temporary workaround to get the NGS data, till we have it in the database.
# To create the RDS files, look into r/scripts/ICWincuNGSdata/incu_incellNGS_GEP00001.R script to recreate them
data <- readRDS("data/GEP00001_data.RDS")
data <- data[data$Allele.fraction > 0.15, c('Plate', 'Well', 'guide', 'Alleles', 'Indel.length',
                                            'Allele.fraction', 'Chromosome', 'Position',
                                            'Variant.type.consequence', 'Symbol..Gene.ID.')]
data$Content <- paste0('MCF7 clone3 ', data$guide)
data <- data[order(data$Plate, data$Well), ]

NGSdata <- readRDS('data/GEP00001_dataNGS.RDS') %>%
   mutate('Content' = as.factor(paste0('MCF7 clone3 ', Layout)),
           'Plate' = as.factor(Plate)) %>%
   mutate('Plate' = gsub('plate', '', Plate)) %>%
   droplevels %>%
   arrange(desc(Plate, Well))

NGSdata.cells <- NGSdata[grepl('-C', NGSdata$Sample),] # extracted cells only
NGSdata.gDNA <- NGSdata[grepl('-G', NGSdata$Sample),]  # extracted gDNA only
#/%DB%

# --------------------------------------------------------------------------------
# connect to database
db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

# connect to tables
cell_line <- tbl(db, "cell_line")        %>% rename(cell_line_id=id, cell_line_name=name, cell_line_description=description)
project <- tbl(db, "project")            %>% rename(project_id=id, project_geid=geid, project_description=description)
experiment_layout <- tbl(db, "experiment_layout") %>% rename(experiment_layout_id=id, experiment_layout_geid=geid)
clone <- tbl(db, "clone")                %>% rename(clone_id=id, clone_name=name)
well <- tbl(db, "well")                  %>% rename(well_id=id)
well_content <- tbl(db, "well_content")  %>% rename(well_content_id=id)
plate <- tbl(db, "plate")                %>% rename(plate_id=id)
growth <- tbl(db, "growth")              %>% rename(growth_id=id)
abundance <- tbl(db, "abundance")        %>% rename(abundance_id=id)
target <- tbl(db, "target")              %>% rename(target_id=id, target_name=name)
guide <- tbl(db, "guide")                %>% rename(guide_id=id, guide_name=name)
guide_well_content_association <- tbl(db, "guide_well_content_association")
# TODO add NGS tables

create_classifier <- function(cellline, clone, guide, content)
{
  thing <- paste(cellline, clone, guide, content, sep=" ")
  thing <- gsub(" NA", "", thing)
  thing <- gsub(" sample", "", thing)
  thing
}

# read protein abundance data
protein_abundance_data <- left_join(well, abundance, by="well_id") %>%
  left_join(., well_content, by="well_content_id") %>%
  left_join(., clone, by="clone_id") %>%
  left_join(., cell_line, by="cell_line_id") %>%
  left_join(., plate, by="plate_id") %>%
  left_join(., guide_well_content_association, by="well_content_id") %>%
  left_join(., guide, by="guide_id") %>%
  left_join(., target, by="target_id") %>%
  collect %>%
  filter(!is.na(content_type)) %>%
  mutate(Well = paste0(row, column)) %>%
  mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
  select(classifier, plate_id, Well, intensity_channel_700, intensity_channel_800) %>%
  rename(Content=classifier, Plate = plate_id, Channel700=intensity_channel_700, Channel800=intensity_channel_800) %>%
  filter(!grepl('background', Content)) %>%            #the background control (no antibodies is not meaningful for proteini amount)
  mutate('ratio800to700' = Channel800/Channel700) %>%  #calculate relative protein abundance, ratio 800 to 700
  mutate(Plate = as.factor(Plate), Well = as.factor(Well))  #change Plate and position to factors, so grouping in ggplot can work

# read clone growth data
clone_growth_data <- left_join(well, growth, by="well_id") %>%
  inner_join(., well_content, by="well_content_id") %>%
  inner_join(., clone, by="clone_id") %>%
  inner_join(., cell_line, by="cell_line_id") %>%
  inner_join(., plate, by="plate_id") %>%
  left_join(., guide_well_content_association, by="well_content_id") %>%
  left_join(., guide, by="guide_id") %>%
  left_join(., target, by="target_id") %>%
  collect %>%
  mutate(position=paste0(row, column)) %>%
  mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
  #select(plate_id, target_name, guide_name, position, hours, confluence_percentage) %>%
  rename(Content = classifier, Plate=plate_id, Target=target_name, Guide=guide_name, Well=position, Elapsed=hours, Confluence=confluence_percentage) %>%
  filter(!is.na(Confluence)) %>%
  mutate(Plate = as.factor(Plate), Well = as.factor(Well), Guide = as.factor(Guide))
#add reactive filters (sliders) to filter out clones with no or little growth

# --------------------------------------------------------------------------------
# processing data

## Convert clone_growth_data to matrix to use with the grofit package
# Trim wells with growth increase (time 93 - time 0) < 10%
growthincrease <- clone_growth_data[clone_growth_data$Elapsed == 93, 'Confluence'] - clone_growth_data[clone_growth_data$Elapsed == 0, 'Confluence']
growthincrease <- data.frame(clone_growth_data[clone_growth_data$Elapsed == 0, c('Plate', 'Well')],
                        growthincrease) 
#it could be worth it to make the growthincrease threshold manually selectable (so a user can say 'I prefer faster growers, etc')
growthincrease <- growthincrease[growthincrease$Confluence > 12,] %>%
mutate(index = paste0(Plate, Well))

clone_growth_data_rate <- select(clone_growth_data, Plate, Well, Content, Elapsed, Confluence) %>%
  mutate(classifier = as.factor(Content)) %>%
  filter(paste0(clone_growth_data$Plate, clone_growth_data$Well) %in% growthincrease$index) %>%
  dcast(Plate + Well + Content ~ Elapsed, value.var = 'Confluence')

clone_growth_data_time <- matrix(data = unique(clone_growth_data$Elapsed), 
                                 ncol = length(unique(clone_growth_data$Elapsed)), 
                                 nrow = nrow(clone_growth_data_rate), byrow = T)

clone_growth_data_grofit <- grofit(clone_growth_data_time, clone_growth_data_rate, control = grofit.control(interactive=F))
clone_growth_data_grofitsum <- summary(clone_growth_data_grofit$gcFit) %>%
  select(1:9, 13)
colnames(clone_growth_data_grofitsum)[c(1:3, 9, 10)] <- c('Plate', 'Well', 'Content', 'mu', 'stdmu') #mu is the maximum slope, and stdmu its std deviation
# NOTE: sometimes the model can't be fitted and we get NA's in mu and stdmu, which are not plotted. 
# Curves look ragged but fine, so we are losing data here due to lack of model fitting

# merge ICW and incu data in a single table for plotting:
clone_growth_protein <- merge(clone_growth_data_grofitsum, protein_abundance_data, all = T, by = c('Plate', 'Well', 'Content')) %>%
  mutate('minusSD' = mu - stdmu, 'plusSD' = mu + stdmu) %>%
  subset(!grepl('normalisation|background', Content))

# --------------------------------------------------------------------------------
# Should be loaded from DB
#%DB%
fun.by <- function(x) {NGS <- by(x[,], INDICES = list(x$Plate, x$Well, x$Type), FUN = as.data.frame)
                       NGS <- NGS[!do.call(c, lapply(NGS, is.null))]
                       NGS <- lapply(NGS, function(a) data.frame(
                                          'Plate' = a$Plate,
                                          'Well' = gsub('P.', '', a$Well),
                                          'Content' = paste0('MCF7 clone3 ', a$Layout),
                                          'Type' = a$Type,  #SNV or INDEL
                                          'NROW' = nrow(a), 
                                          'homo' = a$AlleleFraction >0.85,
                                          'muthet' = a$AlleleFraction > 0.35 & a$AlleleFraction <0.85 & nrow(a) > 1,
                                          'mutwt' = a$AlleleFraction > 0.35 & a$AlleleFraction <0.85 & nrow(a) == 1,
                                          'has.offtargets' = length(unique(a$Gene)) > 1)
                                          )
                       NGS2 <- lapply(NGS, function(a) {  #consolidate rows of the dataframe
                                          a$homo <- all(a$homo)
                                          a$muthet <- as.logical(sum(a$muthet))   #if it has three rows, two T and one F (e.g. it's a muthet with an off target), it will show as muthet
                                          a$mutwt <- as.logical(sum(a$mutwt))
                                          a$has.offtargets <- all(a$has.offtargets)
                                          return(a)
                                })
                       NGS <- lapply(NGS, function(a) a[1,]) #select first row only
}

fun.by2 <- function(x) {NGS <- by(x[,], INDICES = list(x$Plate, x$Well), FUN = as.data.frame)
           NGS <- NGS[!do.call(c, lapply(NGS, is.null))]}

NGSdataby.cells <- lapply(fun.by(NGSdata.cells), function(a) data.frame(a, 'DNAsource' = 'cells'))
NGSdataby.gDNA <- lapply(fun.by(NGSdata.gDNA), function(a) data.frame(a, 'DNAsource' = 'gDNA'))
NGSdataby <- rbind(do.call(rbind, NGSdataby.cells),
                   do.call(rbind, NGSdataby.gDNA))
NGSdataby$Zygosity <- apply(NGSdataby[,c('homo', 'muthet', 'mutwt')], 1, function(a) c('homo', 'muthet', 'mutwt')[a]) %>%
  lapply(function(a) {if(length(a) == 0) a <- NA; return(a)}) %>%
  do.call(rbind, .) %>%
  as.factor

# to check odd results
  # m <- NGSdataby[duplicated(paste0(NGSdataby$Plate, NGSdataby$Well)), c('Plate', 'Well')]
  # m <- paste0(m$Plate, m$Well)
  # NGSdataby$m2 <- as.character(paste0(NGSdataby$Plate, NGSdataby$Well))
  # subset(NGSdataby, m2 %in% m) %>% arrange(Plate, Well) %>% group_by(Plate, Well)
  # subset(NGSdata.cells, Well == 'P3D7')
  # subset(NGSdataby, m2 == '3D7') %>% arrange(Plate, Well) %>% group_by(Plate, Well)

  # NGSdataby.cells[do.call(rbind, lapply(NGSdataby.cells, function(a) nrow(a) > 2))]
  # fun.by2(NGSdata.cells)[do.call(rbind, lapply(NGSdataby.cells, function(a) nrow(a) > 2))]

ICWincuNGS <- merge(clone_growth_protein, NGSdataby, by = c('Plate', 'Well', 'Content'), all = T)
ICWincuNGS$plusSD <- ICWincuNGS$mu + ICWincuNGS$stdmu
ICWincuNGS$minusSD <- ICWincuNGS$mu - ICWincuNGS$stdmu
# remove background and normalisation controls
ICWincuNGS <- subset(ICWincuNGS, !grepl('normalisation|background', Content)) %>% droplevels
#NGSdata <- subset(NGSdata, !grepl('normalisation|background', Content))
#/%DB%

# --------------------------------------------------------------------------------
# plots
# NOTE: plot x axis labels are cut off, this can be workaround by extending the margings 
  #l <- plotly_build(g)
  #l$layout$margin$b <- l$layout$margin$b + 50
##NOTE: hjust does not translate properly into plotly

# plot ratio against cell line
protein_abundance_plot <- ggplot(protein_abundance_data, aes(x = Content, y = ratio800to700, group = Plate:Well)) +
  # geom_violin(aes(group = NULL), color = 'darkgreen') + # it gives error when transformed with ggploly - unknown column 'Content'
  geom_point() + # position = position_dodge()) messes up with hover in ggplotly, so the wrong wells are shown
  theme_bw() +
  ggtitle('Protein content (InCellWestern)') +
  xlab('Cell line') +
  ylab('Relative protein abundance (ratio 800 to 700 nm)') +
  theme(axis.text.x = element_text(angle=15, vjust=0.5, size=10))

# plot clone growth curves
clone_growth_data$cLayout <- as.factor(sapply(clone_growth_data$Content, function(a) strsplit(a, " ")[[1]][3]))
clone_growth_curve <- ggplot(clone_growth_data, 
                             aes(x = Elapsed, y = Confluence,group = cLayout:Plate:Well,
                             color = cLayout)) +
  geom_line() +
  theme_bw() +
  ggtitle('Clone growth curves') +
  xlab('Time (h)') +
  ylab('Confluence (%)') +
  scale_color_discrete(name="Content")

# plot clone growth rates
clone_growth_rate <- ggplot(subset(clone_growth_data_grofitsum, !is.na(mu)), aes(x = Content, y=mu, group = Plate:Well)) +
  # geom_violin(aes(group = NULL)) + # it gives error when transformed with ggploly - unknown column 'Content'
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = mu - stdmu, ymax= mu + stdmu), position = position_dodge(0.6), size = 0.4) +
  theme_bw() +
  ggtitle('Clone growth rate') +
  xlab('Clone') +
  ylab('maximum growth slope') +
  theme(
    # axis.title.x = element_text(face="bold", colour="#990000", size=20),
    axis.text.x  = element_text(angle=15, vjust=0.5, hjust = 0.5, size=10))
# BIG NOTE position = position_dodge()) messes up with hover in ggplotly, so the wrong wells are shown!!

# combined plot (protein + growth slopes + NGS)
ICWincuNGS_plotly <- plot_ly(ICWincuNGS, type = 'scatter', mode = 'markers',
                             # Hover text:
                             text = ~paste("Plate:Well ", Plate, Well)) %>%
  add_trace(x = ~ratio800to700, y = ~mu, symbol = ~Content, color = I("black")) %>% #error bars not showing in the right place... error_y = ~list(type = "data", array = stdmu)
  add_trace(subset(ICWincuNGS, !is.na(Zygosity)), x = ~ratio800to700, y= ~mu, color = ~Zygosity,
            text = ~paste("Plate:Well:Zygosity ", Plate, Well, Zygosity)) %>%
  layout(xaxis = list(title ='protein abundance'), yaxis = list(title = 'maximum growth slope'))

# NGS plots
# proportion of indels and snvs
# initial number of clones screened is hard-coded. Obtain from corresponding number of library samples

NGS.mutations <- plot_ly(NGSdata.cells, type = 'histogram', x = ~Type) %>%
                  layout(xaxis = list(title = 'Type of mutation'), yaxis = list(title = '% of all mutations'))
  
NGS.variants <-  plot_ly(mutate(NGSdata.cells, Variant = c('frameshift', 'synonymous', 'inframe.del', 'intron.noncoding', 'stop.frameshift')[
                                match(NGSdata.cells$Variant, c('frameshift', 'synonymous','inframe_deletion', 'intron,non_coding_transcript', 'stop_gained,frameshift'))]),
                         type = 'histogram', x = ~Variant) %>%
                  layout(xaxis = list(title = 'Type of variant'),
                         yaxis = list(title = '% of all variants'),
                         margin = list(l = 150))

NGS.zygosity <- plot_ly(mutate(NGSdataby, has.offtargets = c('in-target', 'off-target')[match(NGSdataby$has.offtargets, c(F, T))]),
                        type = 'histogram', x = ~Zygosity, color = ~has.offtargets) %>%
                        layout(xaxis = list(title = 'Zygosity'), yaxis = list(title = '%'))

NGS.distancetocutsite <- plot_ly(data, type = 'histogram', x = ~Indel.length) %>%
                  layout(xaxis = list(title = 'Indel length (nt)'), yaxis = list(title = '%'))

# --------------------------------------------------------------------------------
# Notes
# TODO create one plot per well with plotly
# TODO calculate slope of growth curves with grofit
# TODO Moving average: https://en.wikipedia.org/wiki/Moving_average
# TODO Kernel smoother: https://en.wikipedia.org/wiki/Kernel_smoother
