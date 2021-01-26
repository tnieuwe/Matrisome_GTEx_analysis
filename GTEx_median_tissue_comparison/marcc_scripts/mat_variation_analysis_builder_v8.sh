#!/bin/bash
#The third shell step in gtex analysis
#Create marcc master run file


#Run all tissues will invoke each run on its own
echo '#!/bin/sh' > ~/work2/tnieuwe1/projects/matrisome/gtex_v8/variation_analysis/run_all_mat_v8.sh

#Reading both files together

while read -r general && read -r tissue <&4
do
	#generate r script replacing general and tissue with appropriate labels
	sed "s/general/$general/g" ~/work2/tnieuwe1/projects/matrisome/gtex_v8/variation_analysis/general_mat_variation_analysis_v8.R | sed "s/tish/$tissue/g" > ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_analysis_v8.R

	#Create individual run files
			#Using cat to create multiline variables needed for running code
	cat > ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_load_v8.sh <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --mem=30GB
#SBATCH --mail-type=end
#SBATCH --mail-user= tnieuwe1@jhu.edu
#SBATCH --job-name=general_mat_var
#SBATCH --time=20
#SBATCH --account=mhalush1@jhu.edu

module load R
EOF
#Add the line with variable and executed code
	echo "R --no-save < ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_analysis_v8.R" >> ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_load_v8.sh

#Make specific run file for individual troubleshooting
echo '#!/bin/sh' > ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_run_v8.sh
echo "sbatch ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_load_v8.sh" >> ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_run_v8.sh
#Make the run file executable
chmod +x ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_run_v8.sh

#Annotate master run file
echo "sbatch ~/work2/tnieuwe1/projects/gtex_v8/$general/${general}_mat_variation_load_v8.sh" >> ~/work2/tnieuwe1/projects/matrisome/gtex_v8/variation_analysis/run_all_mat_v8.sh
chmod +x  ~/work2/tnieuwe1/projects/matrisome/gtex_v8/variation_analysis/run_all_mat_v8.sh

done < ~/work2/tnieuwe1/data/gtex_v8/gen_tiss_lists/general_list_test.txt 4< ~/work2/tnieuwe1/data/gtex_v8/gen_tiss_lists/tissue_list_test.txt
