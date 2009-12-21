#!/bin/bash

#for table in chromsequence ; do \
# expdp arolfe/tiebreech@psrg DIRECTORY=AROLFE_BACKUP DUMPFILE=${table}.dmp CONTENT=DATA_ONLY TABLES=${table} && \
# impdp core/jeptarict@psrg DIRECTORY=AROLFE_BACKUP DUMPFILE=${table}.dmp REMAP_SCHEMA=arolfe:core CONTENT=DATA_ONLY ; done


#for table in arraydesign galfiles probedesign probetm probelocation fragdist experiment exptToGenome data ipmeta hybmeta scanmeta mleanalysis mleparameters mleanalysisinputs mleToGenome mleresults bayesanalysis bayesparameters bayesanalysisinputs bayesToGenome bayesresults rosettaanalysis rosettaparameters rosettaanalysisinputs rosettaresults rosettaToGenome ; do \
for table in bayesresults rosettaparameters rosettaanalysisinputs rosettaresults  ; do \
impdp ylchipchip/etebrreaji@psrg DIRECTORY=AROLFE_BACKUP DUMPFILE=${table}.dmp REMAP_SCHEMA=arolfe:ylchipchip CONTENT=DATA_ONLY ; done
# expdp arolfe/tiebreech@psrg DIRECTORY=AROLFE_BACKUP DUMPFILE=${table}.dmp CONTENT=DATA_ONLY TABLES=${table} ; done