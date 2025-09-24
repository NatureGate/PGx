

version 1.0
#You need to declaration version information(version 1.0)
workflow PGx{
  input{
    #You need to define"input"code blocks,otherwise it will cause a run error
    File vcf_file
    File bam_file
    File hla_file
    
  }

  call cyriusTask{
    input:
        bam_file = bam_file,
        hla_file = hla_file
  }

  call pharmCatTask{
    input:
        file = vcf_file
  }
  
  #String filename=basename(file,".vcf")
  output{
   #You need to define the workflow output to the “output” code block
    
    File match_report = pharmCatTask.match_report
    File phenotype_result = pharmCatTask.phenotype_result
    File report_file = pharmCatTask.report_file
  }
  
 }
task pharmCatTask{
  input {
    
    File file
    String filename=basename(file,".vcf")
  }
  command {
    echo "hello,world "
    echo ${file},${filename}
    # docker run --rm \
    pharmcat_pipeline --output-dir . ${file}
    ls -l
    pwd
  }
  runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc_external/longrui_4605c9b9e91b4812a4a9f696ec47ca7a_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 2
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "4Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
  output {
    File match_report = "${filename}.match.json"
    File phenotype_result = "${filename}.phenotype.json"
    File report_file = "${filename}.report.html"
    #pharmcat.example.report.html
  }
}
task cyriusTask{
    input {
        File bam_file
        File hla_file
    }
    command {
        echo "hello,world "
        echo ~{bam_file}>bam_path.txt
        # docker run --rm \
        cyrius --manifest bam_path.txt \
                       --genome 38 \
                       --prefix cyp2d6 \
                       --outDir . \
                       --threads 8
        cyp2d6_results_handler.py --input_file cyp2d6.tsv
        hla_results_handler.py --input_file ~{hla_file}
        cat cyp2d6_result.txt>outcall.txt
        echo >> outcall.txt
        cat hla_result.txt>>outcall.txt
    }
    runtime {
        #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
        docker_url: "stereonote_hpc/longrui_4c21ee03ea5d4b2dabd1539a777685ea_private:latest"
        #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
        req_cpu: 8
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "16Gi"
        #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
    }
    output {
        File outcall_file = "outcall.txt"

    }
}