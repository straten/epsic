
set terminal postscript eps enhanced color "times-roman" 20 size 10, 10 
set output "obliquity.eps"

set view 80, 45

splot [-6:6][-6:6] "000.txt" using 2:3:4, "060.txt" using 2:3:4, "120.txt" using 2:3:4, "160.txt" using 2:3:4, "180.txt" using 2:3:4

unset output

