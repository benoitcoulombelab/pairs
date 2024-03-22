nextflow.enable.dsl = 2

params.output = "$baseDir"

log.info """\
         AF2Complex feature detection
         ===================================
         fasta   : ${params.fasta}
         output  : ${params.output}
         account : ${params.account}
         """
        .stripIndent()

workflow {
  fasta = Channel.fromPath(params.fasta, checkIfExists: true)
  feature_detection(fasta)
}

process feature_detection {
  publishDir "${params.output}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(fasta)

  output:
  file("features/*")

  script:
  """
  mkdir -p features
  bash $baseDir/af2c_fea.sh ${fasta} features
  """
}
