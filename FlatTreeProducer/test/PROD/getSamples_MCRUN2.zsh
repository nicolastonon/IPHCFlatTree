#!/bin/env zsh

tag=RunIISpring16MiniAODv2-PUSpring16

wget --no-check-certificate \
--output-document=samples_MCRUN2.txt \
"https://cmsweb.cern.ch/das/request?view=plain&instance=prod%2Fglobal&input=dataset%3D%2F*%2F*${tag}*%2FMINIAODSIM+|+sort+dataset.name"
