#!/bin/bash

# set -x

N_THREADS=20
MIN_MAPQ=20
BINSIZE=10000

parse_line(){
	local line=$1

	sample_path=$(echo $line | cut -d',' -f 1)
	raw_path=$(echo $line | cut -d',' -f 2)
	raw_path=$(realpath ${raw_path})
	genome_assembly=$(echo $line | cut -d',' -f 3)
	genome_sequence=$(echo $line | cut -d',' -f 4)
	genome_sequence=$(realpath ${genome_sequence})
	chromsizes=$(echo $line | cut -d',' -f 5)
	chromsizes=$(realpath ${chromsizes})
	blacklisted_regions=$(echo $line | cut -d',' -f 6)
	if [[ ! -z ${blacklisted_regions} ]]; then
		blacklisted_regions=$(realpath ${blacklisted_regions})
	fi

	echo "*************************************"
	echo "- Sample path: ${sample_path}"
	echo "- Raw path: ${raw_path}"
	echo "- Assembly: ${genome_assembly}"
	echo "- Genome sequence: ${genome_sequence}"
	echo "- Chromosome sizes: ${chromsizes}"
	echo "- Blacklisted regions: ${blacklisted_regions}"
	echo "*************************************"
}

prepare_sample_path(){
	mkdir -p ${sample_path}
	sample_path="$(realpath ${sample_path})"
	log_file=${sample_path}/log.txt
	> ${log_file}
	fastq_dir="${sample_path}/fastq"
	mkdir -p ${fastq_dir}
	for f in `find "${raw_path}" -name "*.fastq*"`
	do
		if [[ ! -e ${fastq_dir}/$(basename ${f}) ]]; then 
			echo "Linking ${f} to ${fastq_dir}" | tee -a ${log_file}
			ln -s ${f} ${fastq_dir}/
		fi
	done
}

get_r1r2_fastq(){
	r1_path=$(ls -l ${fastq_dir}/*_R1_*.fastq* | awk 'NR==1{print $9}')
	r2_path=$(ls -l ${fastq_dir}/*_R2_*.fastq* | awk 'NR==1{print $9}')
}

run_fastqc(){
	qc_path="${fastq_dir}/qc"
	echo "Running FASTQC"
	if [[ ! -d ${qc_path} ]]; then
		mkdir -p ${qc_path}
		fastqc ${r1_path} ${r2_path} -o ${qc_path}
	fi
}


alignment(){
	aligned_path="${sample_path}/aligned"
	mkdir -p ${aligned_path}

	echo "Alignment and sorting of the reads"
	sorted_bam_path="${aligned_path}/$(basename ${sample_path})_sorted.bam"
	if [[ ! -e ${sorted_bam_path} ]]; then
		bam_path="${aligned_path}/$(basename ${sample_path}).bam"
		if [ ! -e ${bam_path} ]; then
			bwa mem -t ${N_THREADS} ${genome_sequence} ${r1_path} ${r2_path} | \
				samtools view -h -b -q ${MIN_MAPQ} -@ ${N_THREADS} \
				> ${bam_path}
		fi
		samtools sort -@ ${N_THREADS} -o ${sorted_bam_path} ${bam_path}
		rm ${bam_path}
	fi
	
	echo "Crearting index"
	sorted_bam_index_path="${sorted_bam_path}.bai"
	if [[ ! -e ${sorted_bam_index_path} ]]; then
		samtools index ${sorted_bam_path} ${sorted_bam_index_path}
	fi

	echo "Calculating statistics of the reads"
	sorted_bam_stats_path="${aligned_path}/$(basename ${sample_path})_sorted.stats"
	if [[ ! -e ${sorted_bam_stats_path} ]]; then
		samtools stats ${sorted_bam_path} > ${sorted_bam_stats_path}
	fi

	echo "Filtering the reads"
	# Taken from https://github.com/PavriLab/repliseq-nf/blob/master/main.nf
	# Exclude:
	# - not primary alignments (0x0100)
	# - unmapped reads (0x004)
	# - mate unmapped (0x008)
	# Inlcude:
	# - paired reads (0x001)
	filtered_bam_path="${aligned_path}/$(basename ${sample_path})_sorted_filtered.bam"
	if [[ ! -e ${filtered_bam_path} ]]; then
		samtools view -@ ${N_THREADS} \
					  -F 0x0100 \
					  -F 0x004 \
					  -F 0x008 \
					  -f 0x001 \
					  -b \
					  -h \
				${sorted_bam_path} \
				> ${filtered_bam_path}
	fi

	filtered_bam_index_path="${filtered_bam_path}.bai"
	if [[ ! -e ${filtered_bam_index_path} ]]; then
		samtools index ${filtered_bam_path} ${filtered_bam_index_path}
	fi
	
	echo "Removing duplicates from the reads"
	dedup_bam_path="${aligned_path}/$(basename ${sample_path})_sorted_dedup.bam"
	if [[ ! -e ${dedup_bam_path} ]]; then
		samtools rmdup ${filtered_bam_path} ${dedup_bam_path}
	fi

	dedup_bam_index_path="${dedup_bam_path}.bai"
	if [[ ! -e ${dedup_bam_index_path} ]]; then
		samtools index ${dedup_bam_path} ${dedup_bam_index_path}
	fi

	dedup_bam_stats_path="${aligned_path}/$(basename ${sample_path})_sorted_dedup.stats"
	if [[ ! -e ${dedup_bam_stats_path} ]]; then
		samtools stats ${dedup_bam_path} > ${dedup_bam_stats_path}
	fi
}

coverage(){
	echo "Computing coverage tracks"
	signals_path="${sample_path}/signals"
	mkdir -p ${signals_path}

	echo "- Creating bins"
	bins_path="${signals_path}/bins_${BINSIZE}.bed"
	if [[ ! -e ${bins_path} ]]; then
		bedtools makewindows -w ${BINSIZE} -g ${chromsizes} > ${bins_path}
	fi

	echo "- Coverage"
	coverage_path="${signals_path}/$(basename ${sample_path})_coverage.bedGraph"
	if [[ ! -e ${coverage_path} ]]; then
		bedtools intersect -c -a ${bins_path} -b ${dedup_bam_path} > ${coverage_path}
	fi	

	echo "- Filtered coverage (without blacklisted regions)"
	allowed_regions_path="${signals_path}/allowed_regions.bed"
	if [[ ! -e ${allowed_regions_path} ]]; then
		bedtools sort -i ${blacklisted_regions} -g ${chromsizes} | \
				bedtools complement -i stdin -g ${chromsizes} \
				> ${allowed_regions_path}
	fi

	allowed_regions_bam_path="${signals_path}/allowed_regions.bam"
	if [[ ! -e ${allowed_regions_bam_path} ]]; then
		samtools view -@ ${N_THREADS} -L ${allowed_regions_path} -b ${dedup_bam_path} > ${allowed_regions_bam_path}
	fi

	coverage_filtered_path="${signals_path}/$(basename ${sample_path})_coverage_filtered.bedGraph"
	if [[ ! -e ${coverage_filtered_path} ]] && [[ ! -z ${blacklisted_regions} ]]; then
		bedtools intersect -c -a ${bins_path} -b ${allowed_regions_bam_path} > ${coverage_filtered_path}
		rm ${allowed_regions_path}
		rm ${allowed_regions_bam_path}
	fi

	coverate_rpkm_path="${signals_path}/$(basename ${sample_path})_coverage_filtered_RPKM.bedGraph"
	if [[ ! -e ${coverate_rpkm_path} ]]; then
		n_reads=$(samtools view -c ${allowed_regions_bam_path})
		scale_factor=$(echo "1000000/${n_reads}" | bc -l)
		scale_path="${signals_path}/scale_factors.txt"
		echo -e "n_reads\t${n_reads}" > ${scale_path}
		echo -e "scale_factor\t${scale_factor}" >> ${scale_path}

		cat ${coverage_filtered_path} | \
			awk -v scale=${scale_factor} '{print $1,$2,$3,$4*1000*(scale/($3-$2)) }' OFS='\t' \
			> $coverate_rpkm_path
	fi
}

############################

CONFIG_FILE=$1

if [[ -z ${CONFIG_FILE} ]]; then
	echo "Usage: $(basename $0) <config_file>"
	exit -1
fi


line_count=0
while read line
do
	# Skipping header
	line_count=$((line_count+1))
	if [[ ${line_count} -eq 1 ]]; then continue; fi
	parse_line ${line}
	prepare_sample_path
	get_r1r2_fastq
	run_fastqc 
	alignment 
	coverage
done < ${CONFIG_FILE}