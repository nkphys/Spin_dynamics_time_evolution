for i in {10..39}
do 
cp input_dynamics.inp input_dynamics_run.inp
sed -i -e "s/VALUESTATE/${i}/g" input_dynamics_run.inp

time ./dynamics 1orb_MCMF input_dynamics_run.inp
echo "${i} microstate done"

done
