#!/usr/bin/env Rscript

#
# Copyright (C) 2017-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

library('DESeq2')

args = commandArgs(trailingOnly = TRUE)
path_count = args[1]
path_condition = args[2]
path_output = args[3]
p_adjust = args[4]

if (length(args) == 5) {
    path_factor = args[5]
}

# Input
countdata = read.csv(path_count, header=TRUE, row.names=1, check.names=FALSE)
coldata = read.csv(path_condition, header=TRUE, row.names=1, check.names=FALSE)
condition = coldata[,1]
if (length(args) == 5) {
    sfactor = read.csv(path_factor, header=TRUE)
}

# Import
dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

if (length(args) == 5) {
    sizeFactors(dds) = as.numeric(as.matrix(sfactor[2,]))
}

# DESeq
dds <- DESeq(dds)

# Get differential expression results
res = results(dds, pAdjustMethod=p_adjust, independentFiltering=FALSE)

# Merge with normalized count data & Write
resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
write.csv(resdata, file=path_output, row.names=FALSE)
