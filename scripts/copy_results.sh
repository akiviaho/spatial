#!/bin/bash

for SAMPLE in BPH_651  BPH_688  CRPC_278  CRPC_489  CRPC_697  PC_03_6712  PC_15420OIK  PC_4980  PC_7875OIK

do
	rsync -av --exclude='analysis' --exclude='metrics_summary.csv' --exclude='*.cloupe' --exclude='*.bam' --exclude='*bam.bai' /lustre/scratch/kiviaho/prostate_spatial/results/${SAMPLE}/outs/* /lustre/scratch/kiviaho/prostate_spatial/single-cell-mapping/data/Visium/${SAMPLE}

done
