#!/usr/bin/env Rscript

#
# Copyright © 2017 Charles E. Vejnar
#
# Use of this source code is governed by an MIT-style license that can be
# found in the LICENSE file or at https://opensource.org/licenses/MIT.
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
