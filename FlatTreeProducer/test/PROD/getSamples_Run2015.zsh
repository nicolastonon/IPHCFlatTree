#!/bin/env zsh

wget --no-check-certificate \
--output-document=samples_Run2015.txt \
"https://cmsweb.cern.ch/das/request?view=plain&limit=1000&instance=prod%2Fglobal&input=dataset%3D%2F*%2F*Run2015B-PromptReco*%2FMINIAOD+|+sort+dataset.name"
