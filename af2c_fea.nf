nextflow.enable.dsl = 2

params.db_preset = "uniprot"
params.feature_mode = "monomer+organism+fullpdb"
now = java.time.LocalDate.now()
date_format = java.time.format.DateTimeFormatter.ISO_LOCAL_DATE
params.max_template_date = date_format.format(now)
params.output = "$baseDir"

log.info """\
         AF2Complex feature detection
         ===================================
         fasta             : ${params.fasta}
         db_preset         : ${params.db_preset}
         feature_mode      : ${params.feature_mode}
         max_template_date : ${params.max_template_date}
         output            : ${params.output}
         account           : ${params.account}
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
  bash $baseDir/af2c_fea.sh ${fasta} features ${params.db_preset} ${params.feature_mode} ${params.max_template_date}
  """
}
