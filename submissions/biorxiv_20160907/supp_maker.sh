cat ../../manuscript/suppmat.tex | \
     sed -E "s/(includegraphics.*\{)/\1..\/..\/manuscript\//" | \
     sed "s/externalcitedocument{pooled_error}/externalcitedocument{final}/" | grep -v "\\\usepackage{natbib}" | \
     egrep -v "(begin|end){color}" | sed "s/\\\newcommand\\\revised\[1\]{\\\textcolor{red}{#1}}/\\\newcommand\\\revised[1]{#1}/" | \
     sed "s/\\\cite{lun2016pooling}/Lun \\\textit{et al.} \\\cite{lun2016pooling}/" |    
     sed "s/\\\cite{kolod2015single}/Kolodziejczyk \\\textit{et al.} \\\cite{kolod2015single}/" |
     sed "s/citep{/cite{/g" \
     > Supp.tex 
pdflatex Supp.tex 
pdflatex Supp.tex

