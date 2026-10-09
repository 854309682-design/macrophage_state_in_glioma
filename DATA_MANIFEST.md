# 数据清单（Raw Data Manifest）

本项目遵循「不把原始大文件提交进 git」的原则：所有原始数据统一保存在仓库外的共享目录
`Datasets/`（默认与本仓库同级；可用环境变量 `GLIOMA_DATA` 指向其他位置），
本文件是该目录的数据**清单（manifest）**。

## 原始数据位置

| 数据集 | 类型 | 本地位置 | 用途 | 状态 |
|---|---|---|---|---|
| GSE103224 (Yuan 2018, *Genome Med*) | scRNA-seq (10x) | `Datasets/GSE103224/` + `Datasets/GSE103224_*.tar` | discovery（本地已有） | ✅ 已核实：无 header 基因×细胞 UMI（60,725 基因 × 23,793 细胞） |
| GSE131928 (Neftel 2019, *Cell*) | scRNA-seq | `Datasets/GSE131928_RAW.tar` + `Datasets/GSE131928/`（解压后） | discovery + 恶性注释参考 | ⚠️ RAW.tar 已下，**待解压** |
| GSE70630 (Tirosh 2016, *Nature*) | scRNA-seq | `Datasets/GSE70630/` + `Datasets/GSE70630_OG_processed_data_v2.txt.gz` | discovery（跨 IDH 泛化） | ✅ |
| GSE89567 (Venteicher 2017, *Science*) | scRNA-seq | `Datasets/GSE89567_IDH_A_processed_data.txt.gz` | IDH-mut 星形细胞瘤微环境 | ✅（IDH_A；重复副本已清理） |
| GSE163120 (Pombo Antunes 2021, *Nat Neurosci*) | scRNA/cite-seq | `Datasets/GSE163120_RAW.tar`（含 human 矩阵）+ 顶层 mouse 文件 | H1 验证（human） | ⚠️ human 矩阵在 tar 内，**待解压** |
| Klemm 2020 (*Cell*) | bulk 分选表达 | `Datasets/klemm_cell_2020/`（BrainTIME counts + clinical） | 髓系标记/生存 signature 参考 | ✅（用户经 brainTIME 获取；无 GEO accession） |
| TCGA-GBM / TCGA-LGG | bulk RNA-seq | `Datasets/TCGA/star_counts/*.tsv`（925 文件）+ GDC 临床 `Datasets/TCGA/TCGA_GBM_LGG_cases_clinical.tsv` + **H3 协变量** `Datasets/TCGA/TCGA_lgggbm_clinical.tsv`（cBioPortal `lgggbm_tcga_pub`） | bulk validation #1 | ✅ 已自动下载（GDC，STAR-Counts）；H3 临床改用 cBioPortal（见下） |
| CGGA mRNAseq_325 / mRNAseq_693 | bulk RNA-seq | `Datasets/CGGA/`（clinical + RSEM + Read_Counts 各一套） | bulk validation #2/#3 | ✅ 已自动下载并核实 |
| 10x GBM Visium (FFPE) | Visium | — | spatial discovery | 🔒 待下载（需 10x 账号） |
| Ravi 2022 (*Cancer Cell*) | Visium + snRNA-seq | — | spatial validation | 🔒 待下载（GEO 未定位到该 series，疑 EGA 受控） |
| GICC 胶质瘤 GWAS meta | GWAS summary stats | — | MR/coloc（可选） | 🔒 待下载（GWAS Catalog 无法按 trait 过滤；疑受控访问/需授权） |

## 参考数据库（NicheNet / SCENIC，仓库外）

| 资源 | 位置 | 大小 | 用途 | 状态 |
|---|---|---|---|---|
| NicheNet 人类 `ligand_target_matrix_nsga2r_final.rds` | `Datasets/nichenet/` | 262 MB | Step 5 NicheNet 配体活性（targets×ligands 矩阵） | ✅ 已就位（自 WoundHealing/refs 复制） |
| NicheNet 人类 `lr_network_human_21122021.rds` | `Datasets/nichenet/` | 20 KB | Step 5 NicheNet 配体-受体网络 | ✅ 已就位 |
| NicheNet `weighted_networks_*`（可选） | — | — | 加权 L-R 网络（配体活性不依赖） | ⬜ 未下载 |
| pySCENIC cisTarget hg38 mc9nr `10kb_up_and_down_tss` feather | `Datasets/scenic_resources/` | 1.25 GB | Step 5 pySCENIC `ctx` 排名库 | ✅ 已获得（2026-09-12，官方 sha1 校验通过）；genes 27,090 / motifs 24,453 |
| 同上 `500bp_up_and_100bp_down_tss` feather | `Datasets/scenic_resources/` | 1.27 GB | pySCENIC `ctx` 排名库（第二个搜索空间） | ✅ 已获得（2026-09-12，sha1 校验通过）；genes 27,015 / motifs 24,453 |
| `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` | `Datasets/scenic_resources/` | 103.6 MB | motif→TF 注释 | ✅ 已获得（163,193 行） |
| `allTFs_hg38.txt` | `Datasets/scenic_resources/` | 11.7 KB | TF 列表 | ✅ 已获得（官方版 1,892 TF，覆盖此前 1,797 TF 的旧副本） |

- 下载脚本：`scripts/05_interaction/step09e_download_scenic_resources.sh`（aria2c 断点续传 + 官方**逐文件 `.sha1sum.txt`** 校验 + ctxcore 完整性校验；注意 aertslab **没有**全局 `sha256sum.txt`）。
- **2026-09-12 网络更新**：`resources.aertslab.org` 在本机**已恢复可达**（此前 TCP 超时不可达）；上述 4 个文件已于当日经该脚本下载完成（2.55 GB）并通过官方 sha1 / ctxcore 校验。
- 说明：`Datasets/scenic_resources/` 与 `Datasets/nichenet/` 均在仓库外共享目录，脚本中以绝对路径常量引用。



- **TCGA-GBM/LGG**：GDC API。**Xena（toil.xenahubs.net）返回 403 AccessDenied**，改用 GDC。
  先试 `POST /data` 打包下载，因服务端组装 925 文件超时（HTTP 500）失败；改用**逐文件并行下载**。
  `curl -P8` 在中国到 GDC 的链路上超时严重（单文件 300s 未完成）；改用 **aria2c**（-j6, resume, retry）
  稳定完成：**925/925 文件，3.7 GB**，无空/半成品文件，表头 `# gene-model: GENCODE v36`。
- **CGGA mRNAseq_325/_693**：`https://cgga.org.cn/download?file=...`（URL 已含 Referer=download.jsp）。
  6 个文件（clinical + RSEM + Read_Counts 各一套），~87 MB。已核实（排除 zip 内 `__MACOSX` 条目）：
  `mRNAseq_325` = 55,524 基因 × **325** 样本；`mRNAseq_693` = 55,524 基因 × **693** 样本。
- **临床**：GDC `cases` 接口；注意 GDC 的 `age_at_diagnosis` 以**天**为单位，且存活/死亡字段分列
  （`demographic.*` 与 `diagnoses.N.*` 数组），分析时需按 case 汇总。
- ⚠️ **H3 临床协变量改用 cBioPortal**（`scripts/07_bulk_validation/step07f_tcga_idh_mgmt.py`）：
  GDC `cases` 接口 (a) 不提供 **IDH/MGMT**，(b) **GBM 缺 `days_to_last_follow_up`**（3/293 非空）→
  存活 GBM 病例被全部丢弃，GBM Cox 沦为 **dead-only**（试跑得 205/207 events 的“all-Dead”假象）。
  故 TCGA H3 的 OS / IDH / MGMT / grade 一律取自 cBioPortal 研究 **`lgggbm_tcga_pub`**
  （patient 级 `OS_MONTHS`/`OS_STATUS`/`AGE`；sample 级 `GRADE`/`IDH_STATUS`/`MGMT_PROMOTER_STATUS`），
  缓存为 `Datasets/TCGA/TCGA_lgggbm_clinical.tsv`（`patientId` = TCGA barcode，与矩阵 `case` 一致）。
  API：`GET /api/studies/lgggbm_tcga_pub/clinical-data?clinicalDataType={PATIENT|SAMPLE}&attributeId=…`。
  CGGA 的 IDH/MGMT 本就在其 clinical zip 内（`step07c` 已在用，无需下载）。

## 目录结构

- 原始数据：`../Datasets/`（共享目录，仓库外；按队列平铺，如 `Datasets/CGGA/`、`Datasets/TCGA/`、`Datasets/IvyGAP/` …）。
- 本清单：`DATA_MANIFEST.md`（仓库内）。
- 脚本中的路径常量：`ROOT`（本仓库根）与 `DATASETS`（`Datasets/`）在各 `scripts/**` 顶部**按脚本自身位置推导**，
  可用环境变量 `GLIOMA_PROJ` / `GLIOMA_DATA` 覆盖 —— 换机器/换目录无需改动任何脚本源代码。

## 注意

- 大文件一律放 `Datasets/`（仓库外）。`GSE131928` 与 `GSE163120` 的 RAW.tar 需解压以获取可用矩阵。
- **TCGA expression 到 clinical 的关联**：`star_counts/<file_id>.tsv` 需按 file_id→case 映射（用 GDC 文件元数据）
  才能在后续与临床表合并；下载时未内嵌该映射，分析 Step 7 需先行补齐。
