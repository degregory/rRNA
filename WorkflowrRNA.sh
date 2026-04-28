#! /bin/bash

echo '===================================================================='
max_children=8

# for file in *_R1_001.fastq.gz
# do
	# Sampid=$(echo $file | rev | cut -d "_" -f 3- | rev)
	# if [[ -f ${Sampid}_R2_001.fastq.gz ]]
	# then
		# bash bbmerge.sh in1=$file in2=${Sampid}_R2_001.fastq.gz  out=$Sampid.merge.fq &>> $Sampid.mergestats.txt
	# else
		# echo "unpaired file ${file}"
	# fi

# done

for file in *.merge.fq
do
	{
	echo 0
	Sampid=$(echo $file | cut -d "." -f 1 )

	if [[ $file == *Fish_16* ]]
	then
	echo $Sampid Fish_16S
	cutadapt -e .3 -g ^GACCCTATGGAGCTTTAGAC -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a CGCTGTTATCCCTADRGTAACT'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	
	elif [[ $file == *Mammal_16* ]]
	then
	echo $Sampid Mammal_16S
	cutadapt -e .3 -g ^CGGTTGGGGTGACCTCGGA -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a GCTGTTATCCCTAGGGTAACT'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	
	elif [[ $file == *Reptile_16* ]]
	then
	echo $Sampid Reptile_16S
	cutadapt -e .3 -g ^AGACNAGAAGACCCTGTG -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a CCTGATCCAACATCGAGG'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	
	elif [[ $file == *MiFishE* ]]
	then
	echo $Sampid MiFishE
	cutadapt -e .3 -g ^GTTGGTAAATCTCGTGCCAGC -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a CATAGTGGGGTATCTAATCCTAGTTTG'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	
	elif [[ $file == *MiFishU* ]]
	then
	echo $Sampid MiFishU
	cutadapt -e .3 -g ^GTCGGTAAAACTCGTGCCAGC -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a CATAGTGGGGTATCTAATCCCAGTTTG'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	
	else
	echo $Sampid 'base'
	cutadapt -e .3 -g ^ACTGGGATTAGATACCCC -o ${Sampid}.cut1.fa $file > ${Sampid}.cut.info
	cutadapt -e .3 -a CTAGAGGAGCCTGTTCTA'$' -o ${Sampid}.cut.fa $Sampid.cut1.fa >> ${Sampid}.cut.info
	fi
	python3 /mnt/c/Weekly_MiSeq/dereprr.py ${Sampid}.cut.fa $Sampid.dereprr.fa 8 &>>  ${Sampid}_derepinfo.txt

	usearch11 -unoise3 $Sampid.dereprr.fa -zotus $Sampid.un.fa -tabbedout $Sampid.unoise3.txt -minsize 1
	python /mnt/c/Weekly_MiSeq/anozotu.py $Sampid
	mafft --auto $Sampid.uc.fa > $Sampid.al.fa 
	} &

	my_pid=$$
	children=$(ps -eo ppid | grep -w $my_pid | wc -w)
	children=$((children-1))
	if [[ $children -ge $max_children ]]
	then
		wait -n
	fi
	echo '||||||||||||||||||||||||||||||||||||||||'

done

wait

python /mnt/c/Weekly_MiSeq/rRNA/1MakeBlastee.py
