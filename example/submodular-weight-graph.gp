set terminal pdfcairo
set out "submodular-weight-graph.pdf"
set xlabel "Submodular weight fraction"
set ylabel "Labeled percent"

set yrange [0:1]
#set xrange [.16:.36]

f(x) = a*x + b
g(x) = c*x + d

fit f(x) "fix.stats" using 3:1 via a,b
fit g(x) "hocr.stats" using 3:1 via c,d

plot "fix.stats" using 3:1 title "FGBZ" with points ps 0.2, \
    f(x) title "FGBZ best fit" with lines, \
    "hocr.stats" using 3:1 title "HOCR" with points ps 0.2, \
    g(x) title "HOCR best fit" with lines
