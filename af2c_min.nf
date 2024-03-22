nextflow.enable.dsl = 2

params.model = "model_"
params.output = "$baseDir"

log.info """\
         AF2Complex relaxation
         ===================================
         targets    : ${params.targets}
         model      : ${params.model}
         output     : ${params.output}
         account    : ${params.account}
         """
        .stripIndent()

workflow {
  targets = Channel.fromPath(params.targets, checkIfExists: true)
  relaxation(targets)
}

process relaxation {
  publishDir "${params.output}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(target)

  output:
  file("structures/*")

  script:
  """
  mkdir -p structures
  cp -r "${params.output}/structures/${target.baseName}" structures
  bash $baseDir/af2c_min.sh ${target} structures structures ${params.model}
  """
}
