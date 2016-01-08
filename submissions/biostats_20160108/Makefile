all: pooled_error.dvi suppmat.pdf

pooled_error.dvi: pooled_error.tex failsim_without_raw.eps failsim_with_raw.eps failsim_with_sum.eps power_with.eps quantile_with.eps power_with_zoom.eps references.bib
	latex pooled_error.tex
	bibtex pooled_error
	latex pooled_error.tex
	latex pooled_error.tex 

failsim_without_raw.eps: ../../failsim/without_raw.eps
	cp $< $@

failsim_with_raw.eps: ../../failsim/with_raw.eps
	cp $< $@

failsim_with_sum.eps: ../../failsim/with_sum.eps
	cp $< $@

power_with.eps: ../../power/with.eps
	cp $< $@

power_with_zoom.eps: ../../power/with_zoom.eps
	cp $< $@

quantile_with.eps: ../../quantile/with.eps
	cp $< $@

suppmat.pdf: ../../manuscript/suppmat.pdf
	cp $< $@

references.bib: ../../manuscript/references.bib
	cp $< $@

