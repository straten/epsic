
set terminal postscript eps enhanced mono "times-roman" 20 
set output "correlation.eps"

set xlabel 'Modulation Index'
set ylabel 'Minimum Correlation Coefficient' rotate by 90
unset key

f(x)=(1/(x*x+1)-1)/(x*x)
plot [0:4] f(x)

unset output

