#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process findTeloReads {
    label "TeloID"
    cpus 1

    input:
        file fastq

    output:
        path "telomeric_read_names.list", emit: telo_read_names
    
    """
    seqkit fq2fa $fastq | NCRF ${params.telomotif} --minlength=${params.minMotifLength} --stats=events > telomeres.ncrf
    cat telomeres.ncrf\
    | grep -B 2 --no-group-sep "^${params.telomotif}+" \
    | grep -A 1 --no-group-sep -P "mRatio=9[0-9]|mRatio=100" | grep -v "^#" | sed 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
    | awk '\$2 - \$5 < 100{print \$1}'\
    > telomeric_read_names.list
    cat telomeres.ncrf\
    | grep -B 2 --no-group-sep "^${params.telomotif}-"\
    | grep -A 1 --no-group-sep -P "mRatio=9[0-9]|mRatio=100" | grep -v "^#" | sed 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
    | awk '\$4 < 100{print \$1}'\
    >> telomeric_read_names.list
    """
}

process pullTeloSeqs {
    label "TeloID"
    cpus 1
   
    input: 
      file fastq
      file telo_read_names
    
    output:
      path "TeloID_telomericReads.fastq", emit: TeloReads

  """
  catfishq $fastq \
  | fgrep -f $telo_read_names -A 3 --no-group-sep /dev/stdin \
  > TeloID_telomericReads.fastq
  """
}

process TeloMap {
    label "TeloID"
    cpus params.threads

    input:
      file ref
      file TeloReads

    output:
      path "telomereReads.bam", emit: Telobam
      path "telomereReads.bam.bai", emit: Telobai

    """
    minimap2 -ax map-ont -t $params.threads --secondary=no $ref $TeloReads \
    | samtools sort -@ 20 > telomereReads.bam
    samtools index telomereReads.bam
    """
}

process DepthCalc {
    label "TeloID"
    cpus params.threads

    input:
      file Telobam
      file Telobai

    output:
      path "telodepth.regions.bed.gz", emit: Depthbed

    """
    mosdepth -t $params.threads -b $params.window_size -n telodepth $Telobam
    """
}

process Plotting {
  label "TeloID"
  cpus 1

  input:
    file telodepth

  output:
    path "telomere_map_plot.$params.file_suffix", emit: TeloPlot

  """
  python $workflow.projectDir/scripts/TeloPLOT.py --depthfile $telodepth --windowsize 1 --output telomere_map_plot.$params.file_suffix
  """
}

process OutputLocation {
  label "TeloID"
  cpus 1

  publishDir "${params.out_dir}", mode: 'copy', pattern: "*"

  input:
    file fname

  output:
    file fname

  """
  echo 'Writing output files'
  """

}

//Workflow entrypoints

workflow {

   main: 
    fastq = file(params.fastq, type: "file")
    ref = file(params.ref, type: "file")
    telolist = findTeloReads(fastq) 
    telomericreads = pullTeloSeqs(fastq, telolist)
    Telobam = TeloMap(ref, telomericreads)
    telodepth = DepthCalc(Telobam)
    finalfig = Plotting(telodepth)
    OutputLocation(finalfig)
} 