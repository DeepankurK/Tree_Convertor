for i in 4 5 6 7 8 9 10
do
	cd dataset"$i"
	cp ../dataset1/runD1.slurm runD"$i".slurm
	sed -i "s/dataset1/dataset"$i"/g" runD"$i".slurm
	sed -i "s/d1/d"$i"/g" runD"$i".slurm
	sbatch runD"$i".slurm
	cd ..
done
