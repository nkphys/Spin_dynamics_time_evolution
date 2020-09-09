#!/bin/bash
#Skw_conf_1.txt
#Skw_using_Crt.txt
infile="Dqw_conv0.01_2states.txt"
outfile="Dqw_conv0.01_2states_Path.txt"
tempfile="temp.txt"
lx=8
ly=8
kx=0
ky=0
index=0
rm ${outfile}
lyby2=$(echo "scale=0 ; ${ly}/2" | bc)
lxby2=$(echo "scale=0 ; ${lx}/2" | bc)


#[ (0,0)--------->(pi,0) ]
kx=0
ky=0
while [ $kx -le ${lxby2} ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"

kx=$[$kx+1]
index=$[$index+1]
done




#( (pi,0)---------> (pi,pi) ]
kx=${lxby2}
ky=1
while [ $ky -le ${lyby2} ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"

ky=$[$ky+1]
index=$[$index+1]
done

#( (pi,pi)---------> (0,0) ]
kx=$[${lxby2}-1]
ky=$[${lyby2}-1]
while [ $ky -ge 0 ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"


ky=$[$ky-1]
kx=$[$kx-1]
index=$[$index+1]
done

#Extra index , "use pm3d corners2color c1" in gnuplot
kx=$[${lxby2}-1]
ky=$[${lyby2}-1]
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"


