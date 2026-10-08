metadata：样本信息表，QC_metadata为QC样本的metadata，里面包含internal-QC，QC1（也就是pooled plasma），Blank是用不到的样本；sample_metadata为study sample的metadata.
LH_ratio_table：HELPs强度表，里面包含3个sheet，3个sheet里的行名列名都是一样的，第一个是L，也就是轻标肽段强度表，第二个是H，也就是重标肽段强度表，第三个是L/H的ratio表，列是样本，里面的样本包含了QC和study sample，可能还有一些用不到的样本，所以需要先剔除metadata的sample_id中不存在的样本列和QC metadata里的Blank.
pg_matrix：蛋白强度定量表，前5列为注释列，后续为样本列，里面的样本包含了QC和study sample，可能还有一些用不到的样本，所以需要先剔除metadata的sample_id中不存在的样本列和QC metadata里的Blank.
sample-QC-InternalQC-HELP_yeast.xlsx：人为加入的两条yeast蛋白 ENO1 酶解的肽段在所有样本中的强度。

## Filtered 97-HELP dataset

54 study samples, 6 Internal-QC, and 15 QC1 pooled plasma samples are retained (75 total). Blank and metadata-unlisted samples have been excluded from L, H, ratio_LH, pg_matrix.tsv, and QC metadata. The original 97 peptides and retained intensity/ratio values are unchanged. See FILTERING_REPORT.md.
