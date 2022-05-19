for i in 2 3 4 5 6 7 8 9 10
do
	cd dataset"$i"
	mkdir SCITE_rslts
	cd SCITE_rslts
	python ~/SiCloneFit/helperFiles/varMat2SCITE.py ../noisy_genotype_dataset"$i".txt > SCITE_noisy_genotype_dataset"$i".csv
	cp ../../dataset1/SCITE_rslts/runSCD1.slurm runSCD"$i".slurm
	sed -i "s/dataset1/dataset"$i"/g" runSCD"$i".slurm
	sed -i "s/d1/d"$i"/g" runSCD"$i".slurm
	sbatch runSCD"$i".slurm
	cd ../..
done
