nextflow.enable.dsl = 2

now = java.time.LocalDate.now()
date_format = java.time.format.DateTimeFormatter.ISO_LOCAL_DATE
params.max_template_date = date_format.format(now)
params.output = "$baseDir"

log.info """\
         AlphaFold 3 data pipeline
         ===================================
         json              : ${params.json}
         database          : ${params.database}
         max_template_date : ${params.max_template_date}
         output            : ${params.output}
         account           : ${params.account}
         """
        .stripIndent()

workflow {
  json = Channel.fromPath(params.json, checkIfExists: true)
  data_pipeline(json)
}

process data_pipeline {
  publishDir "${params.output}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(json)

  output:
  file("data/*")

  script:
  """
  mkdir -p data
  cp $baseDir/run_alphafold.py .
  bash $baseDir/alphafold3_data.sh ${json} data ${params.database} ${params.max_template_date}
  """
}
