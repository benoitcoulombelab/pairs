nextflow.enable.dsl = 2

params.outdir = "$baseDir"
params.step = "all"

log.info """\
         PAIRS
         ===================================
         fasta   : ${params.fasta}
         step    : ${params.step}
         outdir  : ${params.outdir}
         """
        .stripIndent()

workflow {
  fasta = Channel.fromPath(params.fasta, checkIfExists: true)
  alphafold(fasta)
}

process alphafold {
  publishDir "${params.outdir}", mode: "copy"
  stageInMode = "symlink"

  input:
  tuple file(fasta)

  output:
  file("alphafold/*")

  script:
  """
    mkdir alphafold
    bash $baseDir/alphafold.sh ${fasta} alphafold
    """
}
