set terminal pdfcairo
set out "energy-graph.pdf"
set xlabel "Time (seconds)"
set ylabel "Energy"
set xrange [0:60]

plot "medium-test-pc.stats" using 7:9 title "Pairwise cover" with linespoints, \
    "medium-test-pc-grid.stats" using 7:9 title "Pairwise cover (grid)" with linespoints, \
    "medium-test-hocr.stats" using 7:9 title "HOCR" with linespoints, \
    "medium-test-fix.stats" using 7:9 title "FGBZ" with linespoints, \
    "medium-test-grd-heur.stats" using 7:9 title "GRD-Heur" with linespoints
