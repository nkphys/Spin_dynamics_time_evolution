#!/bin/bash
#Skw_conf_1.txt
#Skw_using_Crt.txt
# U-2.0_JClassical0.05/16x16_Tmax1000_dT0.01_RKOrder4_Np192/Temp0.0010/Fourier_w_conv0.001_Microstates_No1/Skw_using_Crt_w_conv0.001.txt
#U-2.0_JClassical0.05/16x16_Tmax1000_dT0.01_RKOrder4_Np192_AnsatzStripeRandom0.05/Temp0.0010/Fourier_w_conv0.01_Microstates_No2/
for Np in 192 #34 32 30 28
do
N_MS=20
#U-2.0_JClassical0.05/16x16_Tmax1000_dT0.01_RKOrder4_Np192_AnsatzStripeRandom0.05/Temp0.0010/Fourier_w_conv0.01/Fw_w_conv0.01_state_1.txt
#U-2.0_JClassical0.05/16x16_Tmax4000_dT0.01_RKOrder4_Np192_AnsatzStripeRandom0.05/Temp0.0010/Fourier_w_conv0.002/Fw_w_conv0.002_state_41.txt
infile="Sqw_check2.txt"
#"Sqw_conv0.0_50states_1stwconv_0.0005_2ndwconv0.03_Tmax8000.txt"
#"U-2.0_JClassical0.05/16x16_Tmax8000_dT0.01_RKOrder4_Np192_AnsatzStripeRandom0.05/Temp0.0010/Fourier_w_conv0.0005_threads4asked_but_notused//Sqw_conv0.0005_50states.txt"
outfile="Sqw_check2_Path2.txt"
#"U-2.0_JClassical0.05/16x16_Tmax8000_dT0.01_RKOrder4_Np192_AnsatzStripeRandom0.05/Temp0.0010/Fourier_w_conv0.0005_threads4asked_but_notused//Sqw_conv0.0005_50states_Path2.txt"
tempfile="temp.txt"
lx=8
ly=8
kx=0
ky=0
index=0
rm ${outfile}
lyby2=$(echo "scale=0 ; ${ly}/2" | bc)
lxby2=$(echo "scale=0 ; ${lx}/2" | bc)
lym1=$(echo "scale=0 ; ${ly}-1" | bc)
lxm1=$(echo "scale=0 ; ${lx}-1" | bc)


#[ (0,0)--------->(\pi,0) ]
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




#( (pi,1)---------> (pi,2*pi-1) ]

kx=${lxby2}
ky=1
while [ $ky -le ${lym1} ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"

ky=$[$ky+1]
index=$[$index+1]
done

#( (pi-1,2pi-1)---------> (0,2pi-1) ]
kx=$[${lxby2}-1]
ky=$[${lym1}]
while [ $kx -ge 0 ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"

awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}


cat ${tempfile} >> ${outfile}

echo " " >> "${outfile}"

kx=$[$kx-1]
index=$[$index+1]
done

#( (0,2pi-2)---------> (0,pi) ]
kx=0
ky=$[${lym1}-1]
while [ $ky -ge ${lyby2} ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"
awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}
cat ${tempfile} >> ${outfile}
echo " " >> "${outfile}"
ky=$[$ky-1]
index=$[$index+1]
done


#((1,pi)--------->(2pi-1,pi) ]
kx=1
ky=${lyby2}
while [ $kx -le ${lxm1} ]
do
rm ${tempfile}
grep "^${kx}   ${ky}  " ${infile} > "${tempfile}"
awk -v var="$index" '{print var"\t"$0}' ${tempfile} > ${tempfile}.bk && mv ${tempfile}.bk ${tempfile}
cat ${tempfile} >> ${outfile}
echo " " >> "${outfile}"
kx=$[$kx+1]
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

echo "${Np} done"
done #Np
