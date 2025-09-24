

version 1.0
#You need to declaration version information(version 1.0)
workflow PrePharmcat{
    input {
        File vcf_file
        File vcf__file_index
        
    }

    call prePharmCatTask{
        input:
            vcf_file = vcf_file,
            vcf__file_index = vcf__file_index
            
    }
  
    output{
    #You need to define the workflow output to the “output” code block
        
        File match_report = pharmCatTask.match_report
        File phenotype_result = pharmCatTask.phenotype_result
        File report_file = pharmCatTask.report_file
        
    }
  
 }



task prePharmCatTask{
  input {
    
    File vcf_file
    File vcf__file_index
    
  }
  command {
    mkdir -p results/pharmcat_ready/
    python3 pharmcat_vcf_preprocessor.py \
      -vcf data/PharmCAT_tutorial_get-rm_wgs_30x_grch38.NA18526.vcf.bgz \
      -refFna reference.fna.bgz \
      -refVcf pharmcat_positions.vcf.bgz \
      -o results/pharmcat_ready/
  }
  runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc_external/longrui_8725a03a53424685888b56480b3f7f08_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 2
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "4Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
  output {
    Array[File] all_files = glob("results/pharmcat_ready/*")
    #pharmcat.example.report.html
  }
}
