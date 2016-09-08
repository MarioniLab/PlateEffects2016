cat ../../manuscript/pooled_error.tex | sed -e '1,/\\end{abstract}/d' | cat preamble.tex - | \
     sed "s/section{/section*{/g" | sed "s/bibliographystyle{biorefs}/bibliographystyle{abbrv}/" | 
     sed -E "s/(includegraphics.*\{)/\1..\/..\/manuscript\//" | sed "s/bibliography{references}/bibliography{..\/..\/manuscript\/references}/" |
     egrep -v "(begin|end){color}" | grep -v "section*{Supplementary Material}" | grep -v "Supplementary Material is available online" |
     sed "s/\\\cite{jaitin2014massively}/Jaitin \\\textit{et al.} \\\cite{jaitin2014massively}/" |
     sed "s/\\\cite{klein2015droplet}/Klein \\\textit{et al.} \\\cite{klein2015droplet}/" |    
     sed "s/\\\cite{kolod2015single}/Kolodziejczyk \\\textit{et al.} \\\cite{kolod2015single}/" |
     sed "s/citep{/cite{/g" \
     > final.tex
pdflatex final.tex
bibtex final
pdflatex final.tex
pdflatex final.tex

