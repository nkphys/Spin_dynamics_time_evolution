
L=6
for w_no in {0..30}
do
awk -v W_No=${w_no} '$5==W_No' Skw_using_Crt_w_conv0.02.txt > Cut_omega${w_no}_temp.txt
awk -v n=${L} '1; NR % n == 0 {print ""}' Cut_omega${w_no}_temp.txt > Cut_omega${w_no}.txt
rm Cut_omega${w_no}_temp.txt
done
