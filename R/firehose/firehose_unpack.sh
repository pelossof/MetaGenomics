#!/bin/bash
day=$1
month=$2
year=$3
s=$4

echo "unpacking firehose stddata__${month}/${day}/${year}/${s}"

for s in stddata__${year}_${month}_${day}/*; do
    study=`echo $s | awk -F '/' '{print $2}'`
    echo ${study}
    
    folder=${s}/${year}${month}${day}
    rnaga=${folder}/*illuminaga*Level_3__RSEM_genes_normalized__data.Level_3.*.tar.gz
    rnahi=${folder}/*illuminahiseq*Level_3__RSEM_genes_normalized__data.Level_3.*.tar.gz
    prot=${folder}/*RPPA_AnnotateWithGene.Level_3.*.tar.gz
    clin=${folder}/*.Clinical_Pick_Tier1.Level_4.*.tar.gz
    clin=${folder}/*.Merge_Clinical.Level_1.*.tar.gz
    maf=${folder}/*.Mutation_Packager_Calls.Level_3.*.tar.gz
    mafraw=${folder}/*.Mutation_Packager_Raw_Calls.Level_3.*.tar.gz
    meth=${folder}/gdac.broadinstitute.org_${study}.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar.gz
    cna=${folder}/gdac.broadinstitute.org_${study}.Merge_cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.Level_3.*.tar.gz
    cnag=${folder}/gdac.broadinstitute.org_${study}.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.*.tar.gz
    countga=${folder}/*illuminaga*RSEM_genes__data.Level_3*.tar.gz
    counthi=${folder}/*illuminahiseq*RSEM_genes__data.Level_3*.tar.gz
    
    if [ -f ${rnaga} ]; then
        echo "${rnaga} exists, unpacking..."
        tar xvzf ${rnaga} -C ${folder}
    fi
    tar xvzf ${rnahi}   -C ${folder}
    tar xvzf ${prot}    -C ${folder}
    tar xvzf ${clin}    -C ${folder}
    tar xvzf ${maf}     -C ${folder}
    tar xvzf ${mafraw}     -C ${folder}
    tar xvzf ${meth}    -C ${folder}

    if [ -f ${countga} ]; then
        echo "${rnaga} exists, unpacking..."
        tar xvzf ${countga} -C ${folder}
    fi
    tar xvzf ${counthi} -C ${folder}

    if [ -f ${cna} ]; then
        echo "${cna} exists, unpacking..."
        tar xvzf ${cna} -C ${folder}
    fi

    if [ -f ${cnag} ]; then
        echo "${cnag} exists, unpacking..."
        tar xvzf ${cnag} -C ${folder}
    fi


    maffolder=${folder}/gdac.broadinstitute.org_${study}.Mutation_Packager_Calls.Level_3.${year}${month}${day}00.0.0
    mafout=${maffolder}/${study}.maf
    echo ${maffolder}
    echo ${mafout}
    echo "processing maf files, and grouping to ${mafout}"
    firehose_mergemaf.sh ${maffolder} ${mafout}
    rm ${maffolder}/TCGA*.txt

    # unpack the raw maf files
    maffolder=${folder}/gdac.broadinstitute.org_${study}.Mutation_Packager_Raw_Calls.Level_3.${year}${month}${day}00.0.0
    mafout=${maffolder}/${study}.maf
    echo ${maffolder}
    echo ${mafout}
    echo "processing raw maf files, and grouping to ${mafout}"
    firehose_mergemaf.sh ${maffolder} ${mafout}
    rm ${maffolder}/TCGA*.txt

    methfolder=${folder}/gdac.broadinstitute.org_${study}.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.${year}${month}${day}00.0.0 
    methfile=${methfolder}/${study}.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt
    echo "processing meth file"
    firehose_meth.sh ${methfile}
    rm $methfile

    countfolder=${folder}/gdac.broadinstitute.org_${study}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.${year}${month}${day}00.0.0
    countfile=${countfolder}/${study}.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt
    echo "processing count file"
    firehose_count.sh $countfile
    rm $countfile

    if [ -f ${countga} ]; then
            countfolderga=${folder}/gdac.broadinstitute.org_${study}.Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.${year}${month}${day}00.0.0
            countfilega=${countfolderga}/${study}.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt
            firehose_count.sh $countfilega
            rm $countfilega
    fi

done

