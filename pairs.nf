nextflow.enable.dsl = 2

params.outdir = "$baseDir"
params.step = "all"

log.info """\
         PAIRS
         ===================================
         fasta   : ${params.fasta}
         step    : ${params.step}
         outdir  : ${params.outdir}
         account : ${params.account}
         """
        .stripIndent()

workflow {
  fasta = Channel.fromPath(params.fasta, checkIfExists: true)
  if (params.step != "alphafold") {
    fasta_prepare = prepare_alphafold(fasta)
  }
  if (params.step == "alphafold") {
    fasta_prepare = fasta.map(fa -> new Tuple(fa, file("${params.outdir}/prepare/${fa.baseName}")))
  }
  if (params.step != "prepare") {
    alphafold(fasta_prepare)
  }
}

process prepare_alphafold {
  publishDir "${params.outdir}", mode: "copy"
  stageInMode = "symlink"

  input:
  file(fasta)

  output:
  tuple file("prepare/${fasta}"), file("prepare/${fasta.baseName}")

  script:
  """
  mkdir -p prepare
  cp ${fasta} prepare
  bash $baseDir/alphafold.sh ${fasta} prepare "prepare"
  """
}

process alphafold {
  publishDir "${params.outdir}", mode: "copy"
  stageInMode = "symlink"

  input:
  tuple file(fasta), file(prepare)

  output:
  file("alphafold/*")

  script:
  """
  mkdir -p alphafold
  cp -r ${prepare} alphafold
  rm -r ${prepare}
  bash $baseDir/alphafold.sh ${fasta} alphafold "alphafold"
  """
}
