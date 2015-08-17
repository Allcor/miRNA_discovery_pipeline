#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mydata = read.table(file=args[1], head=TRUE, sep=";")

library(limma)
pdf(args[2])
vennDiagram(mydata)
dev.off()
