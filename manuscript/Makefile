all: pooled_error.pdf suppmat.pdf

pooled_error.pdf: pooled_error.tex failsim_without_raw.pdf failsim_with_raw.pdf failsim_with_sum.pdf power_with.pdf quantile_with.pdf power_with_zoom.pdf references.bib
	pdflatex pooled_error.tex
	bibtex pooled_error
	pdflatex pooled_error.tex
	pdflatex pooled_error.tex 

failsim_without_raw.pdf: ../failsim/without_raw.eps
	epstopdf $< -o=$@

failsim_with_raw.pdf: ../failsim/with_raw.eps
	epstopdf $< -o=$@

failsim_with_sum.pdf: ../failsim/with_sum.eps
	epstopdf $< -o=$@

power_with.pdf: ../power/with.eps
	epstopdf $< -o=$@

power_with_zoom.pdf: ../power/with_zoom.eps
	epstopdf $< -o=$@

quantile_with.pdf: ../quantile/with.eps
	epstopdf $< -o=$@

suppmat.pdf: suppmat.tex pooled_error.pdf quantile_schematic.pdf realdata_top_edgeR_ola1_sum.pdf realdata_top_edgeR_ola2_sum.pdf
	pdflatex suppmat.tex

quantile_schematic.pdf: ../quantile/sim_schematic.pdf
	cp $< $@

realdata_top_edgeR_ola1_sum.pdf: ../realdata/top_edgeR_ola1_sum.pdf
	cp $< $@

realdata_top_edgeR_ola2_sum.pdf: ../realdata/top_edgeR_ola2_sum.pdf
	cp $< $@
