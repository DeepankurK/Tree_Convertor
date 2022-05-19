for i in 1 2 3 4 5 6 7 8 9 10
do
	cd dataset"$i"
	cd SCITE_rslts
	python ~/SiCloneFit/helperFiles/processSCITEgv.py SCITE_noisy_genotype_dataset"$i"_ml0.gv 100 50 
	cd ../..
done
