set terminal pdfcairo
set out "submodular-weight-graph.pdf"
set xlabel "Submodular weight fraction"
set ylabel "Labeled percent"

plot "fix.stats" using 3:1 title "FGBZ" with points, \
    "hocr.stats" using 3:1 title "HOCR" with points
