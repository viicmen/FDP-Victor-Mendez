#!/bin/bash

.  /etc/profile.d/lmod.sh
.  /etc/profile.d/soft-modules.sh

module load Python/2.7.18-GCCcore-10.2.0
module load Modeller/10.3-Python2
module load MEME/5.1.1-GCCcore-10.2.0-Python-2.7.18
module load EMBOSS/6.6.0-foss-2020b
module load ClustalW2/2.1-foss-2020b
module load Clustal-Omega/1.2.4-foss-2020b
module load BLAST/2.11.0-Linux_x86_64
module load HMMER/3.3.2-foss-2020b
module load x3dna/2.3
module load HBPLUS

'''
For each DBD sequence in the non-redundant dataset, use ModCRE to model the structure and then predict their PWM.
'''

for db in 'Cisbp' 'Jaspar' 'Hocomoco';
do
for file in $db/dbds/* ;
do
python /users/sbi/webservices/modcre/scripts/model_protein.py -i $file --pdb=/users/sbi/webservices/modcre/pdb -v --all --dummy=./dummy -o $db/Modcre/pdbs
done
python /users/sbi/webservices/modcre/scripts/pwm_pbm.py -i $db/Modcre/pdbs -o $db/Modcre/pwms --pdb=/users/sbi/webservices/modcre/pdb --pbm=/users/sbi/webservices/modcre/pbm -v --auto --known --parallel --info=$db/Modcre/pwms/pwm_info.log --dummy=./dummy
done

