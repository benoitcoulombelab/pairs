nextflow.enable.dsl = 2

params.features = "$baseDir/features"
params.preset = "super"
params.models = "model_1_multimer_v3," +
    "model_2_multimer_v3," +
    "model_3_multimer_v3," +
    "model_4_multimer_v3," +
    "model_5_multimer_v3"
params.model_preset = "multimer_np"
params.output = "$baseDir"

log.info """\
         AF2Complex structure prediction
         ===================================
         targets      : ${params.targets}
         features     : ${params.features}
         preset       : ${params.preset}
         models       : ${params.models}
         model_preset : ${params.model_preset}
         output       : ${params.output}
         account      : ${params.account}
         """
        .stripIndent()

workflow {
  targets = Channel.fromPath(params.targets, checkIfExists: true)
  structure_prediction(targets)
}

process structure_prediction {
  publishDir "${params.output}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(target)

  output:
  file("structures/*")

  script:
  """
  mkdir -p structures
  ln -s ${params.features} features
  bash $baseDir/af2c_mod.sh ${target} features structures ${params.preset} ${params.models} ${params.model_preset}
  """
}
