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
      catfishq $fastq \
      | seqkit fq2fa \
      | NCRF $params.telomere_sequence --stats=events --minlength=45 \
      | sed -e 's/[A-Za-z]*=//g' -e 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
      | tr -d "%" \
      | awk '/^#/{p=0}/^#/ && \$4>90{p=1} !/^TTAGGG/ && !/^#/ && p>0 && \$2>4000 && ((/TAACCC/&& \$4 < 30)||(/TTAGGG/&&\$2-\$5 < 30)) {print \$1}' \
      > telomeric_read_names.list 
    """
}

//Workflow entrypoints

workflow {

   main: 
    fastq = file(params.fastq, type: "file")
    ref = file(params.ref, type: "file")
    telolist = findTeloReads(fastq) 
    //telomericreads = pullTeloSeqs(fastq, telolist)
    //Telobam = TeloMap(ref, telomericreads)
    //telodepth = DepthCalc(Telobam)
    //finalfig = Plotting(telodepth)
    //OutputLocation(finalfig)
} 