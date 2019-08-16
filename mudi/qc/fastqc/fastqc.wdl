task run_fastqc{
  File fastq
  String? extra_args
  Int memory
  Int disk_space
  Int num_threads
  Int num_preempt

  command {
    set -euo pipefail

    # Run FastQC
    fastqc \
    ${fastq} \
    --threads ${num_threads} \
    ${default="" extra_args} \
    --outdir '.'

  }


  output {
    File fastqc_html=glob("*_fastqc.html")[0]
    File fastqc_zip=glob("*_fastqc.zip")[0]
  }

  runtime {
    docker: "gcr.io/broad-cga-sanand-gtex/fastqc:latest"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: "${num_threads}"
    preemptible: "${num_preempt}"
  }

  meta {
    author: "Shankara Anand"
  }

}

workflow fastqc{
  call run_fastqc
}
