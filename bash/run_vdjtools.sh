
## Init objects to vdjtools scripts

#############
me=$(whoami)
application_files=/Users/$me/Dropbox/tlgll
files=/Users/$me/Dropbox/tlgll/

vdj=$application_files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$application_files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$application_files/data/selected_tcrb/databases/vdjdb_new


#############

cd $files/tcrb_data/unsorted/tlgll/
t_lgll=$(ls -d "$PWD"/*);

cd $files/tcrb_data/unsorted/emerson/
emerson=$(ls -d "$PWD"/*);

cd $files/tcrb_data/unsorted/kerr/
kerr=$(ls -d "$PWD"/*);

cd $files/tcrb_data/unsorted/savola/
savola=$(ls -d "$PWD"/*);

cd $files/tcrb_data/unsorted/yusko/
yusko=$(ls -d "$PWD"/*);

###############
cd $files
###############
clear

## Preprocess; convert into more human readable format
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $t_lgll $files/tcrb_data/unsorted/t_lgll/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $emerson $files/tcrb_data/unsorted/emerson/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $kerr $files/tcrb_data/unsorted/kerr/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $savola $files/tcrb_data/unsorted/savola/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $yusko $files/tcrb_data/unsorted/yusko/

## Remove unproductive clonotypes
java -Xmx4G -jar $vdj FilterNonFunctional $t_lgll $files/tcrb_data/unsorted/t_lgll/
java -Xmx4G -jar $vdj FilterNonFunctional $emerson $files/tcrb_data/unsorted/emerson/
java -Xmx4G -jar $vdj FilterNonFunctional $kerr $files/tcrb_data/unsorted/kerr/
java -Xmx4G -jar $vdj FilterNonFunctional $savola $files/tcrb_data/unsorted/savola/
java -Xmx4G -jar $vdj FilterNonFunctional $yusko $files/tcrb_data/unsorted/yusko/

## Downsample to 30k
java -Xmx4G -jar $vdj DownSample --size 30000 $t_lgll $files/tcrb_data/subsampled/t_lgll/
java -Xmx4G -jar $vdj DownSample --size 30000 $emerson $files/tcrb_data/subsampled/emerson/
java -Xmx4G -jar $vdj DownSample --size 30000 $kerr $files/tcrb_data/subsampled/kerr/
java -Xmx4G -jar $vdj DownSample --size 30000 $savola $files/tcrb_data/subsampled/savola/
java -Xmx4G -jar $vdj DownSample --size 30000 $yusko $files/tcrb_data/subsampled/yusko/

## CalcDiversityStats
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 30000 $t_lgll $emerson $kerr $savola $yusko $files/tcrb_data/results/diversity/30k_sampled
