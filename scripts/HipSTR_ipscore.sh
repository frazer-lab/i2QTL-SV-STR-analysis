#$ -cwd
#$ -V
#$ -l short
#$ -e error_test_histr
#$ -o out_test_hipstr
#$ -pe smp 2


hostname >& 2

reglist=$1
bamfiles=$2
outdir=$3

input=`head -n $SGE_TASK_ID $reglist | tail -n 1`
name=`basename $input ".bed"`
outname=${outdir}/${name}.vcf.gz
ref_ipscore=/frazer01/publicdata/gatk_bundle_2.8/b37/human_g1k_v37_decoy_Sendai.fasta
ref_hipsci=/frazer01/publicdata/hs37d5/hs37d5.fa

# cmd="/frazer01/home/djakubosky/software/HipSTR/HipSTR --bam-files $bamfiles --fasta /frazer01/publicdata/hs37d5/hs37d5.fa --regions $input --str-vcf $outname"

cmd="/frazer01/home/djakubosky/software/HipSTR/HipSTR --bam-files $bamfiles --fasta $ref_ipscore --regions $input --str-vcf $outname"


echo "Executing: $cmd" >& 2
eval $cmd
date >& 2