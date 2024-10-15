#!/usr/bin/Rscript

# 필요한 패키지 로드
library(Seurat)
library(data.table)  # d 객체를 처리하기 위한 패키지

# 디렉토리 내 모든 `.dem.txt.gz` 파일을 처리하는 함수
convert_files_to_RData <- function(input_dir, output_dir) {

  # 입력 디렉토리 내 모든 .dem.txt.gz 파일 리스트 얻기
  files <- list.files(input_dir, pattern = "*.dem.txt.gz", full.names = TRUE)

  # 파일마다 순회하며 처리
  for (file in files) {

    # 파일명 추출 (확장자 제외)
    file_name <- basename(file)
    file_base <- gsub(".dem.txt.gz", "", file_name)

    # GSM 번호 추출
    gsm_number <- as.numeric(gsub("GSM(\\d+)_.*", "\\1", file_base))
    anno_gsm_number <- gsm_number + 1
    anno_file_name <- paste0("GSM", anno_gsm_number, "_", sub("GSM\\d+_(.*)", "\\1", file_base), ".anno.txt.gz")

    # 짝이 맞는 anno.txt.gz 파일 경로 생성
    anno_file <- file.path(input_dir, anno_file_name)

    # anno 파일이 존재하는지 확인
    if (!file.exists(anno_file)) {
      stop(paste("Annotation file not found for:", file_base))
    }

    # 파일 읽기 (압축 해제 및 읽기)
    data <- read.csv(gzfile(file), header = TRUE, row.names = 1, sep = "\t")

    # Seurat 객체로 변환
    seurat_obj <- CreateSeuratObject(counts = data)

    # d 객체 생성
    d <- as.data.frame(as.matrix(data))  # d 객체를 데이터프레임으로 변환

    # d.stats 생성 및 추가
    d.stats <- data.frame(matrix(ncol = 5, nrow = ncol(d)))  # 열 하나 추가
    colnames(d.stats) <- c("cellID", "umi", "umi.ribo", "umi.mito", "umi.ercc")
    rownames(d.stats) <- colnames(d)
    
    # 셀 ID 추가
    d.stats$cellID <- colnames(d)

    # UMI 통계 계산
    d.stats$umi <- colSums(d)
    d.stats$umi.ribo <- colSums(d[is.element(rownames(d), c("RNA5S", "RNA18S", "RNA28S")), ])
    d.stats$umi.mito <- colSums(d[is.element(rownames(d), c("MTRNR2")) | grepl("^MT-", rownames(d)), ])
    d.stats$umi.ercc <- colSums(d[grepl("^ERCC-", rownames(d)), ])

    # RData 파일로 저장
    save(seurat_obj, d, d.stats, file = file.path(output_dir, paste0(file_base, ".star.expr.RData")))

    # 처리된 파일 메시지 출력
    cat(paste("Processed:", file_name, "\n"))
  }
}

# 함수 실행 (입력 디렉토리와 출력 디렉토리 지정)
input_directory <- "path/to/matrix"
output_directory <- "path/to/RData/output"
convert_files_to_RData(input_directory, output_directory)
