#!/bin/bash

## Args:
# 1: username to db
# 2: password
# 3: the "start" date for these GFF entries (today's date, presumably)
# 4: the GFF filename.

java edu.mit.csail.psrg.analysis.DatabaseClient --insertgff --user $1 --pass $2 --species "Saccharomyces cerevisiae" --version "1-25-05" --date $3 $4 
