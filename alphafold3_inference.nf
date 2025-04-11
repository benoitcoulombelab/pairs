nextflow.enable.dsl = 2

params.models = "$HOME/alphafold3-models"
params.output = "$baseDir"

log.info """\
         AlphaFold 3 model inference
         ===================================
         json       : ${params.json}
         database   : ${params.database}
         models dir : ${params.models}
         output     : ${params.output}
         account    : ${params.account}
         """
        .stripIndent()

workflow {
  json = Channel.fromPath(params.json, checkIfExists: true)
  model_inference(json)
}

process model_inference {
  publishDir "${params.output}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(json)

  output:
  file("structures/*")

  script:
  """
  mkdir -p structures
  cp $baseDir/run_alphafold.py .
  bash $baseDir/alphafold3_inference.sh ${json} structures ${params.database} ${params.models}
  """
}
