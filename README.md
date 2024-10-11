# MitoEdit: targeted mitochondrial base editing

## Introduction
Recently developed DdCBEs and TALEDs, composed of the split bacterial toxin DddAtox derived from Burkholderia cenocepacia, DNA-binding TALE domains, and a UGI, with the latter consisting of a TadA protein domain from Escherichia coli instead of a UGI domain, enable C-to-T and A-G conversions, respectively. 
Currently, targeting a specific site on the mitochondrial genome requires empirically testing different spacing regions, which consist of various TALE sequences, making this a long and time-consuming process.
Here, we introduce MitoEdit, a workflow in Python based on the described rules of DdCBEs and TALEDs, to assist scientists in choosing the optimal window for specifically targeting a base, rather than testing each combination possible.

## Overview

This project provides a set of pipelines for processing mitochondrial DNA (mtDNA) sequences for base editing. The tool allows users to analyze sequences and utilize the TALE-NT tool for designing targeted genome editing approaches.

## Features

- **Multiple Pipelines**: Select from various processing pipelines for different experimental setups.
- **Flexible Input**: Accepts mtDNA sequences and additional bystander information.
- **Integration with TALE-NT**: Runs the TALE-NT tool for designing targeted TALE sequences.
- **Output**: Generates FASTA files, Excel files with windows, and processed TALE outputs.

## Installation

