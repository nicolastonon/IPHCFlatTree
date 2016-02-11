#!/bin/env zsh

tag=RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12

wget --no-check-certificate \
--output-document=samples_MCRUN2.txt \
"https://cmsweb.cern.ch/das/request?view=plain&instance=prod%2Fglobal&input=dataset%3D%2F*%2F*${tag}*%2FMINIAODSIM+|+sort+dataset.name"
