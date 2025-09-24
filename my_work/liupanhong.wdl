version 1.0


import "wdl_task/t_anno_hcwgs.wdl" as Annot
import "wdl_task/t_fetchvcf_hcwgs.wdl" as FetchVCF
import "wdl_task/t_nocall_hcwgs.wdl" as Nocall
import "wdl_task/t_pgx_hcwgs.wdl" as PGX
import "wdl_task/t_prs_hcwgs.wdl" as PRS
import "wdl_task/t_single_disease_hcwgs.wdl" as SD


import "WGS_QC.wdl"
import "ancestry_vcf.wdl"
import "Complex_disease.wdl"
import "Feature.wdl"
import "Haplogroup.wdl"
import "gvcf2vcf.wdl"
import "hla.wdl"

workflow BGE_workflow_hcwgs
{
  input {
    String SampleID

    String vcf
    String vcf_index

    String? bamfile
    String? bamfile_index


    String Outdir
    String sex


    
    Boolean use_hla =true
    Boolean use_pgx_nohla =false

    String sourceD="megabolt"   ### megabolt or gatk



    String chenyi_dir = "/Files/sz_history/chenyi"
    String quning_dir = "/Files/sz_history/quning"
    String ref_dir = "/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom"
    String ancestry_bin = "/Files/sz_history/zhenghaihui/ancestry"
    String feature_bin="/Files/sz_history/zhenghaihui/feature"
    String complex_dir="/Files/sz_history/hexiuju/ref_database"

    String lush_ref_dir="/Files/sz_history/zhenghaihui/lush"

    String HLA_dir = "/Files/sz_history/huangfei/BGE/hcWGS/HLA"
    String HLA_ref_dir="/Files/sz_history/huangfei/BGE/database"


    ## deciper
    # fetchvcf
    String fetchvcf_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/vcf"
    String fetchvcf_db_dir = "/Files/sz_history/yeyongbai/personalgenome/lib_new/base_data"
    # nocall
    String nocall_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/nocall"
    String nocall_db_pots = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/database/nocall"
    # pgx
    String pgx_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/pgx"
    String pgx_db_dir = "/Files/sz_history/yeyongbai/personalgenome/lib_new/base_data"
    # annot
    String annot_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/annot"

    File annot_dir1 = "/Files/sz_history/huangfei/BGE/database/annot/transcript"
    File annot_dir2 = "/Files/sz_history/huangfei/BGE/database/annot/annodb"
    File annot_dir3 = "/Files/sz_history/huangfei/BGE/database/annot/tsv"
    # singleDisease
    String sd_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/singleDisease"
    String sd_db_dir = "/Files/sz_history/yeyongbai/personalgenome/lib_new/base_data"
    # prs
    String prs_bin = "/Files/sz_history/wenfang/Pipeline/WDL/lowpass/scripts/prs2"
    String prs_db_dir = "/Files/sz_history/huangfei/BGE/Decipher/PRS"
    # hla
    #String hla_bin = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/scripts/hla"
    #String hla_db_dir = "/Files/sz_history/wenfang/Pipeline/WDL/hcwgs/database/hla"
    #String ref_dir="/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom"

  }

  ## SampleID = SampleName + 唯一标识符
  String OutdirS = "${Outdir}/${SampleID}"
  String SampleName = SampleID

  String docker_img = "stereonote_hpc_external/wenfang_7367261d33584c889c3ae2a8effbff05_private:latest"

  ##GVCF2VCF
  call gvcf2vcf.gvcf2vcf {
      input:
      infile = vcf,
      infile_index = vcf_index,
      sampleID = SampleID,
      lush_ref_dir = lush_ref_dir

  }

  

  ##HLA
#   if ( use_hla ) {

#     call hla.HLA{
#       input:
#         sample = SampleID,
#         bamFile = select_first([bamfile]),
#         bamFile_index = select_first([bamfile_index]),
#         HLA_dir = HLA_dir,
#         sourceD = sourceD,
#         HLA_ref_dir = HLA_ref_dir

#     }

    ## pgx
    # call PGX.t_pgx as pgx {
    #   input:
    #   sample_id          = SampleName,
    #   gender             = sex,
    #   vcf                = fetchvcf.typots_vcf,
    #   hla_typing_file    = HLA.hla,
    #   bin                = pgx_bin,
    #   db_dir             = pgx_db_dir,
    #   outdir             = "${OutdirS}/PGx",
    #   cpu                = 1,
    #   mem                = 1,
    #   docker_image       = docker_img
    # }
#  }

#  if( use_pgx_nohla ) {
#      call PGX.t_pgx_nohla as pgx2 {
#       input:
#       sample_id          = SampleName,
#       gender             = sex,
#       vcf                = fetchvcf.typots_vcf,
#       bin                = pgx_bin,
#       db_dir             = pgx_db_dir,
#       outdir             = "${OutdirS}/PGx",
#       cpu                = 1,
#       mem                = 1,
#       docker_image       = docker_img
#     }
#  }
 
 
#  ##M-Y

#  call Haplogroup.t_HaploGroup_vcf {  # Changed to match imported name
#       input:
#       sample_id = SampleID,
#       in_file = gvcf2vcf.vcf,
#       in_file_index = gvcf2vcf.vcf_index,
#       gender = sex,
#       ref_dir = ref_dir
#   }




#   ##Ancestry
#   call ancestry_vcf.t_ancestry_int_vcf {
#         input:
#         sample_id=SampleID,
#         in_file=gvcf2vcf.vcf,
#         in_file_index=gvcf2vcf.vcf_index,
#         ancestry_bin = ancestry_bin
#     }


#   ##feature
#   call Feature.t_feature {
#         input:
#         sample_id = SampleID,
#         VCFfile = fetchvcf.typots_vcf,
#         gender = sex,
#         feature_bin = feature_bin,

#     }


#   ##complex_disease
#   call Complex_disease.t_complex_disease {
#         input:

#         sample_id = SampleID,
#         VCFfile_merge = gvcf2vcf.vcf,
#         VCFfile_index_merge = gvcf2vcf.vcf_index,
#         gender = sex,
#         ref_dir = ref_dir,
#         complex_dir = complex_dir
#     }


#   # annotation
#   call Annot.t_anno as annot {
#     input:
#       sample_id    = SampleName,
#       vcf          = gvcf2vcf.vcf,
#       vcf_index    = gvcf2vcf.vcf_index,
#       bin          = annot_bin,
#       ref_dir      = ref_dir,
#       annot_dir1   = annot_dir1,
#       annot_dir2   = annot_dir2,
#       annot_dir3   = annot_dir3,
#       outdir       = "${OutdirS}/Annotation",
#       cpu          = 1,
#       mem          = 1,
#       docker_image = docker_img
#   }

#   # small vcf
#   call FetchVCF.t_fetchvcf as fetchvcf {
#     input:
#       sample_id      = SampleName,
#       avcf           = gvcf2vcf.vcf,
#       avcf_index     = gvcf2vcf.vcf_index,
#       fetchvcf_db_dir= fetchvcf_db_dir,
#       bin            = fetchvcf_bin,
#       outdir         = "${OutdirS}/FetchVCF",
#       cpu            = 1,
#       mem            = 1,
#       docker_image   = docker_img
#   }

#   # nocall
#   call Nocall.t_nocall as noCall {
#     input:
#       sample_id    = SampleName,
#       avcf         = fetchvcf.typots_vcf,
#       bin          = nocall_bin,
#       db_pots_dir  = nocall_db_pots,
#       outdir       = "${OutdirS}/Nocall",
#       cpu          = 1,
#       mem          = 1,
#       docker_image = docker_img
#   }



#   # single disease
#   call SD.t_single_disease as singleDisease {
#     input:
#       sample_id     = SampleName,
#       gender        = sex,
#       anno_file     = annot.anno,
#       nocall_file   = noCall.nocall,
#       bin           = sd_bin,
#       db_dir        = sd_db_dir,
#       outdir        = "${OutdirS}/SingleDisease",
#       cpu           = 1,
#       mem           = 1,
#       docker_image  = docker_img
#   }

#   # prs
#   call PRS.t_prs as prs {
#     input:
#       sample_id    = SampleName,
#       vcf          = fetchvcf.typots_vcf,
#       bin          = prs_bin,
#       db_dir       = prs_db_dir,
#       outdir       = "${OutdirS}/PRS",
#       cpu          = 1,
#       mem          = 1,
#       docker_image = docker_img
#   }


  output {


    File vcf_index1 = gvcf2vcf.vcf_index
    File vcf1 = gvcf2vcf.vcf
    
    # File vcf1 = fetchvcf.typots_vcf
    
    # File? result_feature = t_feature.traits_reps
    # Array[File]? result_haplogroup = t_HaploGroup_vcf.y_m_haplogroup
    # File? result_ancestry = t_ancestry_int_vcf.ancestry_reps
    # File? result_complex = t_complex_disease.complex_reps

    # # nocall
    # File nocall_file = noCall.nocall
    # # pgx
    # File? pgx_file = pgx.pgx_reps
    # File? pgx_nocal = pgx.pgx_reps_nocal
    
    # File? pgx_file2 = pgx2.pgx_reps
    # File? pgx_nocal2 = pgx2.pgx_reps_nocal
    # # annotation
    # File anno_file = annot.anno
    # File anno_index_file = annot.anno_index
    # # single disease
    # File? single_disease_file = singleDisease.single_disease
    # # prs
    # File? vd_prs = prs.vitamin_d_prs_reps
    # File? creativity_prs = prs.creativity_prs_reps
    # File? wellbeing_prs = prs.wellbeing_prs_reps
    # File? adhd_prs = prs.adhd_prs_reps
    # File? stroke_prs = prs.stroke_prs_reps
    # # hla
    # File? hla_type = HLA.hla
  }
}