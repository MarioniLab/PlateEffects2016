echo '
source("ESpresso.R")
source("estimateParam.R")
' | Rdevel --no-save

# Upload script.
# rsync -azv results_ESpresso/ cruk:lustre/PlateEffects/reference/results_ESpresso/ --delete
