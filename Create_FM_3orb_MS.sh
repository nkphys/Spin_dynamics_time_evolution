
file_out="MS_FM.txt"
echo "#x              y         n_spin     Theta(x,y)       Phi(x,y)Moment_Size(x,y)" > $file_out

for x in {0..23}
do

for y in {0..23}
do


for spin in 0 1 2
do

a=$((1 + $RANDOM%1000))
a=$(echo "${a}/50000" | bc -l)

echo "${x}     ${y}      ${spin}     ${a}     0       1" >> ${file_out}
done

done

done
