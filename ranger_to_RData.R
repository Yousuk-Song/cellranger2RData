#!/usr/bin/Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

# 첫 번째 인자는 Cell Ranger 출력 디렉토리, 두 번째 인자는 저장할 파일 이름, 없으면 기본값 사용
if (length(args) < 2) {
    stop("Usage: ./ranger_to_Rdata.R <cell_ranger_output_dir> <output_file_name>")
}

cell_ranger_output_dir <- args[1]
output_name <- args[2]

setwd(cell_ranger_output_dir)
generateDigitalExprMatrix <- function(cell_ranger_output_dir, output_name, min.umis=NULL, min.genes=NULL) {
    message("Reading Cell Ranger output from: ", cell_ranger_output_dir)

    seurat_obj <- Read10X(data.dir = cell_ranger_output_dir)
    seurat_obj <- CreateSeuratObject(counts = seurat_obj, min.cells = min.genes, min.features = min.umis)

    if (!is.null(min.umis)) {
        seurat_obj <- subset(seurat_obj, subset = nCount_RNA > min.umis)
    }

    if (!is.null(min.genes)) {
        seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min.genes)
    }

    # Seurat 5.0.0 이후로는 `slot` 대신 `layer`를 사용해야 함
    d <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

    # 데이터프레임으로 변환하기 전에 행과 열 수를 확인
    message("Number of cells: ", length(colnames(d)))
    message("Number of genes: ", nrow(d))

    # UMI, ribosomal, mitochondrial, ERCC stats 계산
    d.stats <- data.frame(cellID = colnames(d), umi = colSums(d))
    d.stats$umi.ribo <- colSums(d[is.element(rownames(d), c("RNA5S", "RNA18S", "RNA28S")), ])
    d.stats$umi.mito <- colSums(d[is.element(rownames(d), c("MTRNR2")) | grepl("^MT-", rownames(d)), ])
    d.stats$umi.ercc <- colSums(d[grepl("^ERCC-", rownames(d)), ])

    # 주어진 경로에 파일 저장
    save(list = c("d", "d.stats"), file = paste0(output_name, ".expr.RData"))

    message("### ", Sys.time(), " - all done\n")
    return(list(data = d, stats = d.stats))
}
# 기본값: 최소 UMI 수와 최소 유전자 수 설정
min.umis <- 1000
min.genes <- 500

# 함수 호출
D <- generateDigitalExprMatrix(cell_ranger_output_dir, output_name, min.umis, min.genes)

