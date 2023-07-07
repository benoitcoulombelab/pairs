nextflow.enable.dsl = 2

params.outdir = "$baseDir"
params.step = "all"

log.info """\
         AlphaFold-pairs
         ===================================
         baits   : ${params.baits}
         targets : ${params.targets}
         step    : ${params.step}
         outdir  : ${params.outdir}
         """
        .stripIndent()

workflow {
  baits = Channel.fromPath(params.baits, checkIfExists: params.step != "alphafold")
  targets = Channel.fromPath(params.targets, checkIfExists: params.step != "alphafold")
  fasta = fasta_pairs(baits, targets) | flatten
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

process fasta_pairs {
  input:
  file(baits)
  file(targets)

  output:
  file("proteins/*.fasta")

  """
  mkdir proteins
  fasta-pairs --baits ${baits} --targets ${targets} --output proteins -u -i
  """
}

process prepare_alphafold {
  publishDir "${params.outdir}", mode: "copy"
  stageInMode = "copy"

  input:
  file(fasta)

  output:
  tuple file("prepare/${fasta}"), file("prepare/${fasta.baseName}")

  script:
  """
  mkdir prepare
  cp ${fasta} prepare
  bash $baseDir/alphafold.sh ${fasta} prepare "prepare"
  """
}

process alphafold {
  publishDir "${params.outdir}", mode: "copy"
  stageInMode = "copy"

  input:
  tuple file(fasta), file(prepare)

  output:
  file("alphafold/*")

  script:
  """
    echo ${fasta}
    echo ${prepare}
    mkdir alphafold
    cp -r ${prepare} alphafold
    rm -r ${prepare}
    bash $baseDir/alphafold.sh ${fasta} alphafold "alphafold"
    """
}
