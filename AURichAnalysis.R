#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# R Data Analysis:
# Pisapia et al (2018)
# Tristetraprolin regulates the turnover of autoimmune-associated HLA-DQ mRNAs 
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/darogan/2018_Pisapia_DelPozzo
#
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics, 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+-------------------------------------------------------------------------------")
message("+ Load the required libraries")
message("+-------------------------------------------------------------------------------")

library("ggplot2")
library("ggrepel")
library("cowplot")
library("reshape")


message("+-------------------------------------------------------------------------------")
message("+ Set up the data directory")
message("+-------------------------------------------------------------------------------")

(DataDir <- getwd())
if (!is.null(DataDir)) setwd(DataDir)
Project    <- "HLA"


message("+-------------------------------------------------------------------------------")
message("+ Create Figure 2A")
message("+-------------------------------------------------------------------------------")


message("+ Figure 2A: Dinucleotide Frequencies")

difreqs          <- read.table("dinuc.table.txt", header=TRUE, stringsAsFactors=TRUE)
difreqs$Sequence <- gsub('DQ', '3DQ', difreqs$Sequence)
difreqs$facet    <- factor(difreqs$Sequence, levels = c("3DQA101", "3DQB105", "3DQA105",  "3DQB102"))

diDQA101_x    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQA101")$Exp_freq)
diDQA105_x    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQA105")$Exp_freq)
diDQA101_y    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQA101")$Obs_freq)
diDQA105_y    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQA105")$Obs_freq)

diDQB102_x    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQB102")$Exp_freq)
diDQB105_x    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQB105")$Exp_freq)
diDQB102_y    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQB102")$Obs_freq)
diDQB105_y    <- mean(subset(difreqs, grepl("UA|AU|UU", Dinucleotide) & Sequence == "3DQB105")$Obs_freq)

multx <- 1
multy <- 1
di.annot      <- data.frame(x=c(diDQA101_x*multx, diDQA105_x*multx, diDQB102_x*multx, diDQB105_x*multx ), 
                            y=c(diDQA101_y*multy, diDQA105_y*multy, diDQB102_y*multy, diDQB105_y*multy ), 
                            facet=c("3DQA101", "3DQA105", "3DQB102", "3DQB105"), 
                            label=c( round(diDQA101_y/diDQA101_x, 2), round(diDQA105_y/diDQA105_x, 2), 
                                     round(diDQB102_y/diDQB102_x, 2), round(diDQB105_y/diDQB105_x, 2)))

plot.di <- ggplot(data=difreqs, aes(x=Exp_freq, y=Obs_freq, colour=Dinucleotide)) + 
  geom_point(alpha=0.25, size=1, colour="grey") +
  geom_point(  data = subset(difreqs, grepl("UA|AU|UU", Dinucleotide)), alpha=0.99, size=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") +
  #geom_text_repel(data = subset(difreqs, grepl("UA|AU|UU", Dinucleotide)), aes(label=Dinucleotide), size=2, nudge_x = 0.05, nudge_y = 0.01) +
  geom_text(data=di.annot, aes(x=0.55,y=0.1,label=paste0(label, " (Exp/Obs)")), inherit.aes=FALSE, size=3) +
  #scale_y_continuous(breaks=seq(0,0.3,0.1)) +
  xlim(0,0.7) + ylim(0,0.3) +
  facet_wrap(~ facet, ncol=2) +
  #xlab("Expected Frequency") + 
  xlab("") + ylab("Observed Frequency") + 
  ggtitle(paste0("Di-nucleotide")) +
  theme(legend.position="none", text = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), 
        strip.text.x=element_text(size=9, colour="black", face="bold"), strip.background=element_rect(colour="white",fill="white")   )

message("+ Figure 2A: Trinucleotide Frequencies")

trifreqs          <- read.table("trinuc.table.txt", header=TRUE, stringsAsFactors=TRUE)
trifreqs$Sequence <- gsub('DQ', '3DQ', trifreqs$Sequence)
trifreqs$facet    <- factor(trifreqs$Sequence, levels = c("3DQA101", "3DQB105", "3DQA105",  "3DQB102"))

triDQA101_x    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQA101")$Exp_freq)
triDQA105_x    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQA105")$Exp_freq)
triDQA101_y    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQA101")$Obs_freq)
triDQA105_y    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQA105")$Obs_freq)

triDQB102_x    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQB102")$Exp_freq)
triDQB105_x    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQB105")$Exp_freq)
triDQB102_y    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQB102")$Obs_freq)
triDQB105_y    <- mean(subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide) & Sequence == "3DQB105")$Obs_freq)
     
multx <- 1
multy <- 1
tri.annot       <- data.frame(x=c(triDQA101_x*multx, triDQA105_x*multx, triDQB102_x*multx, triDQB105_x*multx), 
                              y=c(triDQA101_y*multy, triDQA105_y*multy, triDQB102_y*multy, triDQB105_y*multy), 
                             facet=c("3DQA101", "3DQA105", "3DQB102", "3DQB105"), 
                             label=c( round(triDQA101_y/triDQA101_x, 2), round(triDQA105_y/triDQA105_x, 2), 
                                      round(triDQB102_y/triDQB102_x, 2), round(triDQB105_y/triDQB105_x, 2)))

plot.tri <- ggplot(data=trifreqs, aes(x=Exp_freq, y=Obs_freq, colour=Trinucleotide)) + 
  geom_point(alpha=0.25, size=1, colour="grey") +
  geom_point(  data = subset(trifreqs, grepl("[AU][AU][AU]", Trinucleotide)), alpha=0.99, size=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") +
  geom_text(data=tri.annot, aes(x=0.5,y=0.125,label=paste0(label, " (Exp/Obs)")), inherit.aes=FALSE, size=3) +
  xlim(0,0.6) + ylim(0,0.15) + 
  facet_wrap(~ facet, ncol=2) +
  xlab("") + ylab("Observed Frequency") +
  ggtitle(paste0("Tri-nucleotide")) +
  theme(legend.position="none", text = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), 
        strip.text.x=element_text(size=9, colour="black", face="bold"), strip.background=element_rect(colour="white",fill="white")   )

message("+ Figure 2A: Quadnucleotide Frequencies")

quadfreqs          <- read.table("quadnuc.table.txt", header=TRUE, stringsAsFactors=TRUE)
quadfreqs$Sequence <- gsub('DQ', '3DQ', quadfreqs$Sequence)
quadfreqs$facet    <- factor(quadfreqs$Sequence, levels = c("3DQA101", "3DQB105", "3DQA105",  "3DQB102"))

quadDQA101_x    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQA101")$Exp_freq)
quadDQA105_x    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQA105")$Exp_freq)
quadDQA101_y    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQA101")$Obs_freq)
quadDQA105_y    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQA105")$Obs_freq)

quadDQB102_x    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQB102")$Exp_freq)
quadDQB105_x    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQB105")$Exp_freq)
quadDQB102_y    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQB102")$Obs_freq)
quadDQB105_y    <- mean(subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide) & Sequence == "3DQB105")$Obs_freq)

multx <- 1
multy <- 5
quad.annot      <- data.frame(x=c(quadDQA101_x*multx, quadDQA105_x*multx, quadDQB102_x*multx, quadDQB105_x*multx), 
                              y=c(quadDQA101_y*multy, quadDQA105_y*multy, quadDQB102_y*multy, quadDQB105_y*multy), 
                              facet=c("3DQA101", "3DQA105", "3DQB102", "3DQB105"), 
                              label=c( round(quadDQA101_y/quadDQA101_x, 2), round(quadDQA105_y/quadDQA105_x, 2), 
                                       round(quadDQB102_y/quadDQB102_x, 2), round(quadDQB105_y/quadDQB105_x, 2)))

plot.quad <- ggplot(data=quadfreqs, aes(x=Exp_freq, y=Obs_freq, colour=Quadnucleotide)) + 
  geom_point(alpha=0.25, size=1, colour="grey") +
  geom_point(  data = subset(quadfreqs, grepl("[AU][AU][AU][AU]", Quadnucleotide)), alpha=0.99, size=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") +
  geom_text(data=quad.annot, aes(x=0.4,y=0.05,label=paste0(label, " (Exp/Obs)")), inherit.aes=FALSE, size=3) +
  xlim(0,0.5) + ylim(0,0.06) + 
  facet_wrap(~ facet, ncol=2) +
  xlab("") + 
  ylab("Observed Frequency") +
  ggtitle(paste0("Tetra-nucleotide")) +
  theme(legend.position="none", text = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), 
        strip.text.x=element_text(size=9, colour="black", face="bold"), strip.background=element_rect(colour="white",fill="white")   )


message("+ Figure 2A: Pentanucleotide Frequencies")

pentafreqs          <- read.table("pentanuc.table.txt", header=TRUE, stringsAsFactors=TRUE)
pentafreqs$Sequence <- gsub('DQ', '3DQ', pentafreqs$Sequence)
pentafreqs$facet    <- factor(pentafreqs$Sequence, levels = c("3DQA101", "3DQB105", "3DQA105",  "3DQB102"))

pentaDQA101_x    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQA101")$Exp_freq)
pentaDQA105_x    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQA105")$Exp_freq)
pentaDQA101_y    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQA101")$Obs_freq)
pentaDQA105_y    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQA105")$Obs_freq)

pentaDQB102_x    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQB102")$Exp_freq)
pentaDQB105_x    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQB105")$Exp_freq)
pentaDQB102_y    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQB102")$Obs_freq)
pentaDQB105_y    <- mean(subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide) & Sequence == "3DQB105")$Obs_freq)

multx <- 1
multy <- 5
penta.annot      <- data.frame(x=c(pentaDQA101_x*multx, pentaDQA105_x*multx, pentaDQB102_x*multx, pentaDQB105_x*multx), 
                               y=c(pentaDQA101_y*multy, pentaDQA105_y*multy, pentaDQB102_y*multy, pentaDQB105_y*multy), 
                              facet=c("3DQA101", "3DQA105", "3DQB102", "3DQB105"), 
                              label=c( round(pentaDQA101_y/pentaDQA101_x, 2), round(pentaDQA105_y/pentaDQA105_x, 2), 
                                       round(pentaDQB102_y/pentaDQB102_x, 2), round(pentaDQB105_y/pentaDQB105_x, 2)))

plot.penta <- ggplot(data=pentafreqs, aes(x=Exp_freq, y=Obs_freq, colour=Pentanucleotide)) + 
  geom_point(alpha=0.25, size=1, colour="grey") +
  geom_point(  data = subset(pentafreqs, grepl("[AU][AU][AU][AU][AU]", Pentanucleotide)), alpha=0.99, size=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") +
  geom_text(data=penta.annot, aes(x=0.33,y=0.025,label=paste0(label, " (Exp/Obs)")), inherit.aes=FALSE, size=3) +
  xlim(0,0.4) + ylim(0,0.03) + 
  facet_wrap(~ facet, ncol=2) +
  xlab("") + 
  ylab("Observed Frequency") +
  ggtitle(paste0("Penta-nucleotide")) +
  theme(legend.position="none", text = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), 
        strip.text.x=element_text(size=9, colour="black", face="bold"), strip.background=element_rect(colour="white",fill="white")   )


message("+ Figure 2A: Hexanucleotide Frequencies")

hexafreqs          <- read.table("hexanuc.table.txt", header=TRUE, stringsAsFactors=TRUE)
hexafreqs$Sequence <- gsub('DQ', '3DQ', hexafreqs$Sequence)
hexafreqs$facet    <- factor(hexafreqs$Sequence, levels = c("3DQA101", "3DQB105", "3DQA105",  "3DQB102"))
hexaDQA101_x    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQA101")$Exp_freq)
hexaDQA105_x    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQA105")$Exp_freq)
hexaDQA101_y    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQA101")$Obs_freq)
hexaDQA105_y    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQA105")$Obs_freq)

hexaDQB102_x    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQB102")$Exp_freq)
hexaDQB105_x    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQB105")$Exp_freq)
hexaDQB102_y    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQB102")$Obs_freq)
hexaDQB105_y    <- mean(subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide) & Sequence == "3DQB105")$Obs_freq)

multx <- 1
multy <- 5
hexa.annot      <- data.frame(x=c(hexaDQA101_x*multx, hexaDQA105_x*multx, hexaDQB102_x*multx, hexaDQB105_x*multx), 
                               y=c(hexaDQA101_y*multy, hexaDQA105_y*multy, hexaDQB102_y*multy, hexaDQB105_y*multy), 
                               facet=c("3DQA101", "3DQA105", "3DQB102", "3DQB105"), 
                               label=c( round(hexaDQA101_y/hexaDQA101_x, 2), round(hexaDQA105_y/hexaDQA105_x, 2), 
                                        round(hexaDQB102_y/hexaDQB102_x, 2), round(hexaDQB105_y/hexaDQB105_x, 2)))

plot.hexa <- ggplot(data=hexafreqs, aes(x=Exp_freq, y=Obs_freq, colour=Hexanucleotide)) + 
  geom_point(alpha=0.25, size=1, colour="grey") +
  geom_point(  data = subset(hexafreqs, grepl("[AU][AU][AU][AU][AU][AU]", Hexanucleotide)), alpha=0.99, size=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour="grey") +
  geom_text(data=hexa.annot, aes(x=0.25,y=0.025,label=paste0(label, " (Exp/Obs)")), inherit.aes=FALSE, size=3) +
  xlim(0,0.3) + ylim(0,0.03) + 
  facet_wrap(~ facet, ncol=2) +
  xlab("Expected Frequency") + ylab("Observed Frequency") +
  ggtitle(paste0("Hexa-nucleotide")) +
  theme(legend.position="none", text = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), 
        strip.text.x=element_text(size=9, colour="black", face="bold"), strip.background=element_rect(colour="white",fill="white")   )



message("+-------------------------------------------------------------------------------")
message("+ Use cowlot to make a seperate panel for Figure 2A only ")
message("+-------------------------------------------------------------------------------")

dev.off()
theme_set(theme_cowplot())
pdf(paste0(Project, "_multinuc.pdf") ,width=7,height=13, onefile=FALSE)
par(bg=NA)
plot_grid(plot.di, plot.tri, plot.quad, plot.penta, plot.hexa, labels=c("A", "", "", "", ""), ncol = 1, nrow = 5)
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ Create Figure 2B")
message("+-------------------------------------------------------------------------------")

DQB_densities   <- read.table("DQB105_DQB102_densities.txt", header=TRUE, fill=TRUE, sep=",")
DQB_densities.m <- melt(DQB_densities, id=c("position"))


DQA_densities   <- read.table("DQA101_DQA105_densities.txt", header=TRUE, fill=TRUE, sep=",")
DQA_densities.m <- melt(DQA_densities, id=c("position"))


plot_DQA_di <- ggplot() + 
               geom_jitter(data=DQA_densities, aes(x=position, y=DQA101_di), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQA_densities$DQA101_di > DQA_densities$DQA105_di), "red", "grey")  )  +
               geom_jitter(data=DQA_densities, aes(x=position, y=DQA105_di), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQA_densities$DQA105_di > DQA_densities$DQA101_di), "blue", NA)  )  +
               ylab("DQA*") + xlab("") +
               ylim(0.0001,3.5)

plot_DQB_di <- ggplot() + 
               geom_jitter(data=DQB_densities, aes(x=position, y=DQB105_di), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQB_densities$DQB105_di > DQB_densities$DQB102_di), "red", "grey")  )  +
               geom_jitter(data=DQB_densities, aes(x=position, y=DQB102_di), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQB_densities$DQB102_di > DQB_densities$DQB105_di), "blue", NA)  )  +
               ylab("DQB*") + xlab("") +
               ylim(0.0001,3.5)

plot_DQA_tri <- ggplot() + 
               geom_jitter(data=DQA_densities, aes(x=position, y=DQA101_tri), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQA_densities$DQA101_tri > DQA_densities$DQA105_tri), "red", "grey")  )  +
               geom_jitter(data=DQA_densities, aes(x=position, y=DQA105_tri), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
               col = ifelse( (DQA_densities$DQA105_tri > DQA_densities$DQA101_tri), "blue", NA)  )  +
               ylab("DQA*") + xlab("") +
               ylim(0.0001,3.5)

plot_DQB_tri <- ggplot() + 
                geom_jitter(data=DQB_densities, aes(x=position, y=DQB105_tri), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                col = ifelse( (DQB_densities$DQB105_tri > DQB_densities$DQB102_tri), "red", "grey")  )  +
                geom_jitter(data=DQB_densities, aes(x=position, y=DQB102_tri), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                col = ifelse( (DQB_densities$DQB102_tri > DQB_densities$DQB105_tri), "blue", NA)  )  +
                ylab("DQB*") + xlab("") +
                ylim(0.0001,3.5)

plot_DQA_quad <- ggplot() + 
                 geom_jitter(data=DQA_densities, aes(x=position, y=DQA101_quad), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQA_densities$DQA101_quad > DQA_densities$DQA105_quad), "red", "grey")  )  +
                 geom_jitter(data=DQA_densities, aes(x=position, y=DQA105_quad), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQA_densities$DQA105_quad > DQA_densities$DQA101_quad), "blue", NA)  )  +
                 ylab("DQA*") + xlab("") +
                 ylim(0.0001,3.5)

plot_DQB_quad <- ggplot() + 
                 geom_jitter(data=DQB_densities, aes(x=position, y=DQB105_quad), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQB_densities$DQB105_quad > DQB_densities$DQB102_quad), "red", "grey")  )  +
                 geom_jitter(data=DQB_densities, aes(x=position, y=DQB102_quad), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQB_densities$DQB102_quad > DQB_densities$DQB105_quad), "blue", NA)  )  +
                 ylab("DQB*") + xlab("") +
                 ylim(0.0001,3.5)


plot_DQA_penta <- ggplot() + 
                  geom_jitter(data=DQA_densities, aes(x=position, y=DQA101_penta), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                  col = ifelse( (DQA_densities$DQA101_penta > DQA_densities$DQA105_penta), "red", "grey")  )  +
                  geom_jitter(data=DQA_densities, aes(x=position, y=DQA105_penta), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                  col = ifelse( (DQA_densities$DQA105_penta > DQA_densities$DQA101_penta), "blue", NA)  )  +
                  ylab("DQA*") + xlab("") +
                  ylim(0.0001,3.5)

plot_DQB_penta <- ggplot() + 
                  geom_jitter(data=DQB_densities, aes(x=position, y=DQB105_penta), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                  col = ifelse( (DQB_densities$DQB105_penta > DQB_densities$DQB102_penta), "red", "grey")  )  +
                  geom_jitter(data=DQB_densities, aes(x=position, y=DQB102_penta), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                  col = ifelse( (DQB_densities$DQB102_penta > DQB_densities$DQB105_penta), "blue", NA)  )  +
                  ylab("DQB*") + xlab("") +
                  ylim(0.0001,3.5)


plot_DQA_hexa <- ggplot() + 
                 geom_jitter(data=DQA_densities, aes(x=position, y=DQA101_hexa), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQA_densities$DQA101_hexa > DQA_densities$DQA105_hexa), "red", "grey")  )  +
                 geom_jitter(data=DQA_densities, aes(x=position, y=DQA105_hexa), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQA_densities$DQA105_hexa > DQA_densities$DQA101_hexa), "blue", NA)  )  +
                 ylab("DQA*") + xlab("Sequence Position") +
                 ylim(0.0001,3.5)

plot_DQB_hexa <- ggplot() + 
                 geom_jitter(data=DQB_densities, aes(x=position, y=DQB105_hexa), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQB_densities$DQB105_hexa > DQB_densities$DQB102_hexa), "red", "grey")  )  +
                 geom_jitter(data=DQB_densities, aes(x=position, y=DQB102_hexa), size = 3, shape=15, alpha=0.5, width = 0.025, height=0.01, 
                 col = ifelse( (DQB_densities$DQB102_hexa > DQB_densities$DQB105_hexa), "blue", NA)  )  +
                 ylab("DQB*") + xlab("") +
                 ylim(0.0001,3.5)

title.di    <- ggdraw() + draw_label("Di-nucleotide motif position",    fontface='bold')
title.tri   <- ggdraw() + draw_label("Tri-nucleotide motif position",   fontface='bold')
title.quad  <- ggdraw() + draw_label("Tetra-nucleotide motif position", fontface='bold')
title.penta <- ggdraw() + draw_label("Penta-nucleotide motif position", fontface='bold')
title.hexa  <- ggdraw() + draw_label("Hexa-nucleotide motif position",  fontface='bold')
xaxis.lab   <- ggdraw() + draw_label("Sequence Position", size=10)


grid_di    <- plot_grid( plot_DQA_di,    plot_DQB_di,    labels=c(""), ncol = 2, nrow = 1)
grid_tri   <- plot_grid( plot_DQA_tri,   plot_DQB_tri,   labels=c(""), ncol = 2, nrow = 1)
grid_quad  <- plot_grid( plot_DQA_quad,  plot_DQB_quad,  labels=c(""), ncol = 2, nrow = 1)
grid_penta <- plot_grid( plot_DQA_penta, plot_DQB_penta, labels=c(""), ncol = 2, nrow = 1)
grid_hexa  <- plot_grid( plot_DQA_hexa,  plot_DQB_hexa,  labels=c(""), ncol = 2, nrow = 1)

message("+-------------------------------------------------------------------------------")
message("+ Use cowlot to make a seperate panel for Figure 2B only ")
message("+-------------------------------------------------------------------------------")

dev.off()
pdf(paste0(Project, "_multimotif.pdf") ,width=7,height=13, onefile=FALSE)
par(bg=NA)
plot_grid(title.di, grid_di, title.tri, grid_tri, title.quad, grid_quad, title.penta, grid_penta, title.hexa, grid_hexa,
          ncol = 1, nrow = 11, labels=c("B", "", "", "", "", "", "", "", "", ""),
          rel_heights = c(.1, .5, .1, .5, .1, .5, .1, .5, .1, .5))
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ Use cowlot to make the complete Figure 2")
message("+-------------------------------------------------------------------------------")

ti_grid_di    <- plot_grid( title.di,    grid_di,    labels=c(""), ncol = 1, nrow = 2, rel_heights = c(.1, .9) )
ti_grid_tri   <- plot_grid( title.tri,   grid_tri,   labels=c(""), ncol = 1, nrow = 2, rel_heights = c(.1, .9) )
ti_grid_quad  <- plot_grid( title.quad,  grid_quad,  labels=c(""), ncol = 1, nrow = 2, rel_heights = c(.1, .9) )
ti_grid_penta <- plot_grid( title.penta, grid_penta, labels=c(""), ncol = 1, nrow = 2, rel_heights = c(.1, .9) )
ti_grid_hexa  <- plot_grid( title.hexa,  grid_hexa,  labels=c(""), ncol = 1, nrow = 2, rel_heights = c(.1, .9) )


dev.off()
pdf(paste0(Project, "_multimotif_combo.pdf") ,width=18,height=15, onefile=FALSE)
par(bg=NA)
plot_grid(plot.di,    ti_grid_di, 
          plot.tri,   ti_grid_tri,
          plot.quad,  ti_grid_quad,
          plot.penta, ti_grid_penta,
          plot.hexa,  ti_grid_hexa,
          ncol = 2, nrow = 5, labels=c("A", "B", "", "", "", "", "", "", "", "", ""), label_size = 20 )

dev.off()

message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")
