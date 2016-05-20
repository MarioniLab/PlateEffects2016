all: pooled_error.dvi suppmat.pdf

pooled_error.dvi: pooled_error.tex main_failsim.eps main_okaysim.eps main_power.eps main_power_zoom.eps main_swapped.eps references.bib
	latex pooled_error.tex
	bibtex pooled_error
	latex pooled_error.tex
	latex pooled_error.tex 

pooled_error.tex: ../../manuscript/pooled_error.tex
	cat $< | sed "s/\.pdf\}/.eps}/" > $@

references.bib: ../../manuscript/references.bib
	cp $< $@

suppmat.pdf: ../../manuscript/suppmat.pdf
	cp $< $@

main_failsim.eps: ../../simulations/ESpresso/results_failsim/main_raw.eps
	cp $< $@
	
main_okaysim.eps: ../../simulations/ESpresso/results_failsim/main_sum.eps
	cp $< $@

main_power.eps: ../../simulations/ESpresso/results_power/ROC_1.eps
	cp $< $@

main_power_zoom.eps: ../../simulations/ESpresso/results_power/ROC_zoom_1.eps
	cp $< $@

main_swapped.eps: ../../realdata/all_scrambled.eps
	cp $< $@

