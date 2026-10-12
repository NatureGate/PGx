#!/usr/bin/env python3
"""
从 bcftools 提取的 sequenced VCF 更新 PharmCAT 模板基因型。

规则：
  1. 只保留测序位点 (FORMAT 含 DP, INFO 无 RAF=)，忽略 impute 位点
  2. DP >= 10: 使用 GT;  DP < 10: 不更新, 保持 0/0
  3. 假阳性过滤: 测序位点 DP 达标后, 检查 INFO 注释
     (QD/MQ/FS/SOR/MQRankSum/ReadPosRankSum) + 杂合 AB,
     触发任一条件 → 视为假阳性, 不更新, 模板保持 0/0
  4. 按 CHROM+POS+REF+ALT 匹配 PharmCAT 模板 (多等位基因任一匹配即可)
  5. 匹配且样本有变异 (非 0/0): 更新模板基因型
  6. 其他位点保持不变, 所有基因型均为 0/0 或实际变异

输入:
  --seq        测序VCF (bcftools 提取结果)
  --template   PharmCAT模板VCF
  --sample     样本名称
  --outdir     输出目录

输出:
  - {sample}.pharmcat.seq_updated.vcf       — 更新后的 PharmCAT VCF
  - {sample}.pharmcat.seq_changed.tsv       — 变更记录 (简化格式)
  - {sample}.pharmcat.seq_fp_filtered.vcf   — FP过滤后的干净测序VCF
  - {sample}.outcall.tsv                    — MT-RNR1 变异报告
  - {sample}.annotated.tsv                  — 变异位点注释 (同 annotated.csv 格式)
  - {sample}.nocall.tsv                     — NoCall 位点报告

用法:
  python pgx_seq_update.py --seq seq.vcf --template template.vcf --sample SAMPLE --outdir /data/out
  python pgx_seq_update.py --seq seq.vcf --template template.vcf --sample SAMPLE --outdir /data/out --dp 8
  python pgx_seq_update.py --seq seq.vcf --template template.vcf --sample SAMPLE --outdir /data/out --no-fp-filter
"""

import sys
import os
import re
import time

# ======================== 配置 ========================
SEQ_VCF = None
TEMPLATE_VCF = None
PCR_TARGET = None          # PCR target BED 文件
DP_THRESHOLD = 10  # DP >= 10 才更新

# 假阳性过滤 (仅测序位点, 触发任一条件 -> 视为假阳性, 不更新)
FP_FILTER_ENABLED = True
FP_QUAL_MIN = 20.0       # QUAL < 20.0 → FP
FP_QD_MIN = 3.0          # QD < 2.0  → FP
FP_MQ_MIN = 40.0         # MQ < 40.0 → FP
FP_FS_MAX = 60.0         # FS > 60.0 → FP
FP_SOR_MAX = 3.0         # SOR > 3.0 → FP
FP_MQ_RANK_MIN = -12.5   # MQRankSum < -12.5 → FP
FP_READ_POS_RANK_MIN = -8.0  # ReadPosRankSum < -8.0 → FP
FP_AB_MIN = 0.2          # 杂合 AB < 0.2 → FP
FP_AB_MAX = 0.8          # 杂合 AB > 0.8 → FP

# NoCall 检测
NOCALL_QUAL_MIN = 20.0   # QUAL < 20 → NoCall

outdir = None

# 从文件名推断样本名
SAMPLE_NAME = None
MTRNR1_VAR=[]
MTRNR1_GENE='MT-RNR1'


# 命令行参数解析
USAGE = """用法:
  python pgx_seq_update.py --seq <sequenced_vcf> --template <template_vcf> --pcr-target <bed> --sample <sample_name> --outdir <output_dir>

必填参数:
  --seq        测序VCF文件路径 (bcftools提取结果)
  --template   PharmCAT模板VCF文件路径
  --pcr-target PCR target BED文件路径
  --sample     样本名称
  --outdir     输出目录

可选参数:
  --dp         DP阈值 (默认: 10)
  --no-fp-filter  关闭假阳性过滤
"""

args = sys.argv[1:]
i = 0
while i < len(args):
    if args[i] == "--seq" and i + 1 < len(args):
        SEQ_VCF = args[i + 1]; i += 2
    elif args[i] == "--template" and i + 1 < len(args):
        TEMPLATE_VCF = args[i + 1]; i += 2
    elif args[i] == "--pcr-target" and i + 1 < len(args):
        PCR_TARGET = args[i + 1]; i += 2
    elif args[i] == "--sample" and i + 1 < len(args):
        SAMPLE_NAME = args[i + 1]; i += 2
    elif args[i] == "--outdir" and i + 1 < len(args):
        outdir = args[i + 1]; i += 2
    elif args[i] == "--dp" and i + 1 < len(args):
        DP_THRESHOLD = int(args[i + 1]); i += 2
    elif args[i] == "--no-fp-filter":
        FP_FILTER_ENABLED = False; i += 1
    elif args[i] in ("-h", "--help"):
        print(USAGE)
        sys.exit(0)
    else:
        i += 1

# 检查必填参数, 缺失则交互输入
if SEQ_VCF is None:
    SEQ_VCF = input("请输入测序VCF文件路径 (--seq): ").strip()
if TEMPLATE_VCF is None:
    TEMPLATE_VCF = input("请输入PharmCAT模板VCF文件路径 (--template): ").strip()
if PCR_TARGET is None:
    PCR_TARGET = input("请输入PCR target BED文件路径 (--pcr-target): ").strip()
if SAMPLE_NAME is None:
    SAMPLE_NAME = input("请输入样本名称 (--sample): ").strip()
if outdir is None:
    outdir = input("请输入输出目录 (--outdir): ").strip()

# 校验文件存在
for label, path in [("测序VCF", SEQ_VCF), ("模板VCF", TEMPLATE_VCF), ("PCR target", PCR_TARGET)]:
    if not os.path.isfile(path):
        print(f"错误: {label}文件不存在: {path}")
        sys.exit(1)

# 创建输出目录
os.makedirs(outdir, exist_ok=True)

OUT_VCF = os.path.join(outdir, f"{SAMPLE_NAME}.pharmcat.seq_updated.vcf")
OUT_CHANGED = os.path.join(outdir, f"{SAMPLE_NAME}.pharmcat.seq_changed.tsv")
OUT_FP_FILTERED = os.path.join(outdir, f"{SAMPLE_NAME}.pharmcat.seq_fp_filtered.vcf")
OUTCALL_TSV=os.path.join(outdir, f"{SAMPLE_NAME}.outcall.tsv")
OUT_ANNOTATED = os.path.join(outdir, f"{SAMPLE_NAME}.annotated.tsv")
OUT_NOCALL = os.path.join(outdir, f"{SAMPLE_NAME}.nocall.tsv")

# ======================== 假阳性检测 ========================
def _parse_info(info_str):
    """解析 INFO 字段为 dict"""
    d = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def _get_float(info_dict, key):
    val = info_dict.get(key)
    if val is None or val is True:
        return None
    try:
        return float(val)
    except ValueError:
        return None


def check_false_positive(info_str, fmt_str, sample_str, gt, qual_str):
    """
    检查一个测序位点是否为假阳性。
    返回 (is_fp, reasons)
    """
    if not FP_FILTER_ENABLED:
        return False, []

    reasons = []
    info = _parse_info(info_str)

    # QUAL 过滤
    try:
        qual = float(qual_str) if qual_str != "." else 0.0
    except ValueError:
        qual = 0.0
    if qual < FP_QUAL_MIN:
        reasons.append(f"QUAL={qual:.1f}<{FP_QUAL_MIN:.0f}")

    # INFO 级别过滤
    qd = _get_float(info, "QD")
    mq = _get_float(info, "MQ")
    fs = _get_float(info, "FS")
    sor = _get_float(info, "SOR")
    mq_rank = _get_float(info, "MQRankSum")
    rp_rank = _get_float(info, "ReadPosRankSum")

    if qd is not None and qd < FP_QD_MIN:
        reasons.append(f"QD={qd:.2f}<{FP_QD_MIN}")
    # if mq is not None and mq < FP_MQ_MIN:
    #     reasons.append(f"MQ={mq:.2f}<{FP_MQ_MIN}")
    # if fs is not None and fs > FP_FS_MAX:
    #     reasons.append(f"FS={fs:.2f}>{FP_FS_MAX}")
    # if sor is not None and sor > FP_SOR_MAX:
    #     reasons.append(f"SOR={sor:.2f}>{FP_SOR_MAX}")
    # if mq_rank is not None and mq_rank < FP_MQ_RANK_MIN:
    #     reasons.append(f"MQRankSum={mq_rank:.2f}<{FP_MQ_RANK_MIN}")
    # if rp_rank is not None and rp_rank < FP_READ_POS_RANK_MIN:
    #     reasons.append(f"ReadPosRankSum={rp_rank:.2f}<{FP_READ_POS_RANK_MIN}")

    # 等位基因平衡 (仅杂合)
    fmt_fields = fmt_str.split(":")
    sample_fields = sample_str.split(":")
    fmt_idx = {f: j for j, f in enumerate(fmt_fields)}

    if "GT" in fmt_idx and "AD" in fmt_idx:
        alleles = re.split(r"[/|]", gt)
        unique = set(a for a in alleles if a != ".")
        if len(unique) > 1 and "0" in unique:
            try:
                ad_vals = sample_fields[fmt_idx["AD"]].split(",")
                ad_ref = int(ad_vals[0])
                ad_alt = sum(int(v) for v in ad_vals[1:])
                total = ad_ref + ad_alt
                if total > 0:
                    ab = ad_alt / total
                    if ab < FP_AB_MIN or ab > FP_AB_MAX:
                        reasons.append(f"AB={ab:.3f} outside [{FP_AB_MIN},{FP_AB_MAX}]")
            except (ValueError, IndexError):
                pass

    return len(reasons) > 0, reasons


# ======================== Step 1: 读取 sequenced VCF ========================
def read_sequenced_vcf():
    """
    读取 sequenced VCF, 只保留测序位点 (非 impute)。
    返回两个 dict:
      variant_data: (chrom, pos, ref, alt) → (gt, dp, info, fmt, sample, qual)
      ref_data:     (chrom, pos) → (dp, info, fmt, sample, qual)
    """
    print("=" * 60)
    print("Step 1: 读取 sequenced VCF")
    print("=" * 60)
    print(f"  文件: {SEQ_VCF}")

    variant_data = {}  # (chrom, pos, ref, alt) → (gt, dp, info, fmt, sample, qual)
    ref_data = {}      # (chrom, pos) → (dp, info, fmt, sample, qual)
    n_total = 0
    n_sequenced = 0
    n_imputed = 0
    n_dp_low = 0
    n_ref_calls = 0
    n_fp_filtered = 0
    n_written = 0
    fp_details = []    # 前 20 条假阳性详情

    # 打开 FP 过滤后的输出 VCF
    fout_fp = open(OUT_FP_FILTERED, "w", encoding="utf-8")

    with open(SEQ_VCF, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                # 写 header 到 FP 过滤输出
                fout_fp.write(line)
                continue
            n_total += 1
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            info = fields[7]
            fmt = fields[8]
            sample = fields[9] if len(fields) > 9 else ""

            # 判断是否为 impute 位点: INFO 含 RAF= 或 FORMAT 含 DS
            is_imputed = ("RAF=" in info) or ("DS" in fmt.split(":"))
            if is_imputed:
                n_imputed += 1
                continue

            # 测序位点
            n_sequenced += 1
            fmt_fields = fmt.split(":")
            sample_fields = sample.split(":")
            fmt_idx = {f: j for j, f in enumerate(fmt_fields)}

            # 提取 DP
            dp = 0
            if "DP" in fmt_idx:
                try:
                    dp = int(sample_fields[fmt_idx["DP"]])
                except (ValueError, IndexError):
                    dp = 0

            # 提取 GT
            gt = "0/0"
            if "GT" in fmt_idx:
                try:
                    gt = sample_fields[fmt_idx["GT"]]
                except IndexError:
                    gt = "0/0"

            # DP < threshold → 跳过, 不存储
            if dp < DP_THRESHOLD:
                n_dp_low += 1
                continue

            # 假阳性过滤
            is_fp, fp_reasons = check_false_positive(info, fmt, sample, gt, qual)
            if is_fp:
                n_fp_filtered += 1
                if len(fp_details) < 20:
                    fp_details.append((chrom, pos, ref, alt, gt, dp, fp_reasons))
                continue

            # 通过过滤, 写入 FP 过滤输出
            fout_fp.write(line)
            n_written += 1

            # ALT="." 表示参考纯合 (无变异), 按 (chrom, pos) 存储
            if alt == ".":
                n_ref_calls += 1
                key = (chrom, pos)
                if key not in ref_data or dp > ref_data[key][0]:
                    ref_data[key] = (dp, info, fmt, sample, qual)
            else:
                # 有具体 ALT, 多等位基因 (逗号分隔) 拆成单独条目
                for single_alt in alt.split(","):
                    single_alt = single_alt.strip()
                    if not single_alt:
                        continue
                    key = (chrom, pos, ref, single_alt)
                    if key not in variant_data or dp > variant_data[key][1]:
                        variant_data[key] = (gt, dp, info, fmt, sample, qual)

    fout_fp.close()

    print(f"  总数据行:       {n_total}")
    print(f"  测序位点:       {n_sequenced}")
    print(f"  Impute 位点:    {n_imputed} (已忽略)")
    print(f"  DP < {DP_THRESHOLD}:      {n_dp_low} (跳过)")
    if FP_FILTER_ENABLED:
        print(f"  假阳性过滤:     {n_fp_filtered} (GT→不更新, 保持0/0)")
        print(f"  通过过滤写入:   {n_written}")
        if fp_details:
            print(f"\n  假阳性详情 (前 {len(fp_details)} 条):")
            for c, p, r, a, g, d, reasons in fp_details:
                print(f"    {c}:{p} {r}>{a} GT={g} DP={d} | {'; '.join(reasons)}")
            if n_fp_filtered > len(fp_details):
                print(f"    ... 共 {n_fp_filtered} 条")
            print()
    print(f"  参考纯合 (ALT='.'): {n_ref_calls}")
    print(f"  变异位点 (有ALT): {len(variant_data)}")
    print(f"  唯一参考位点:   {len(ref_data)}")
    print(f"  FP过滤后VCF:    {OUT_FP_FILTERED}")
    print()
    return variant_data, ref_data


# ======================== Step 2: 更新模板 ========================
def update_template(variant_data, ref_data):
    """
    读取 PharmCAT 模板, 用 sequenced 数据更新基因型。
    匹配逻辑:
      1. 先按 (chrom, pos, ref, alt) 匹配变异位点
      2. 再按 (chrom, pos) 匹配参考纯合位点 (ALT=".")
    """
    print("=" * 60)
    print("Step 2: 更新 PharmCAT 模板")
    print("=" * 60)
    print(f"  模板: {TEMPLATE_VCF}")
    print(f"  输出: {OUT_VCF}")
    print()

    n_total = 0
    n_matched_variant = 0
    n_matched_ref = 0
    n_updated = 0
    n_unchanged = 0
    n_not_found = 0
    changed_records = []

    with open(TEMPLATE_VCF, "r", encoding="utf-8") as fin, \
         open(OUT_VCF, "w", encoding="utf-8") as fout, \
         open(OUT_CHANGED, "w", encoding="utf-8") as fchg, \
         open(OUT_ANNOTATED, "w", encoding="utf-8") as fann:

        fchg.write("#CHROM\tPOS\tREF\tALT\tID\tGENE\tOLD_GT\tNEW_GT\tDP\n")
        fann.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSM\n")

        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            rsid = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            old_gt = fields[9] if len(fields) > 9 else "0/0"

            # 提取基因名
            gene = "."
            for kv in info.split(";"):
                if kv.startswith("PX="):
                    gene = kv[3:]
                    break

            n_total += 1
            matched = False
            gt_updated = False  # 标记本次是否发生了基因型更新

            # 优先匹配变异位点 — 模板 ALT 可能是多等位基因 (逗号分隔)
            # 只要测序中任一 ALT 匹配模板中任一 ALT 即可
            seq_matched_data = None  # 存储匹配到的测序数据 (info, fmt, sample, qual)
            for single_alt in alt.split(","):
                single_alt = single_alt.strip()
                if not single_alt:
                    continue
                var_key = (chrom, pos, ref, single_alt)
                if var_key in variant_data:
                    new_gt, dp, seq_info, seq_fmt, seq_sample, seq_qual = variant_data[var_key]
                    seq_matched_data = (seq_info, seq_fmt, seq_sample, seq_qual)
                    n_matched_variant += 1
                    matched = True
                    if new_gt != old_gt and new_gt != "0/0":
                        fields[9] = new_gt
                        gt_updated = True
                        if gene==MTRNR1_GENE:
                            MTRNR1_VAR.append(f'm.{pos}{ref}>{alt}')
                            
                        changed_records.append((chrom, pos, ref, alt, rsid, gene, old_gt, new_gt, dp))
                        fchg.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rsid}\t{gene}\t{old_gt}\t{new_gt}\t{dp}\n")
                        n_updated += 1
                    else:
                        n_unchanged += 1
                    break

            # 再匹配参考纯合位点 (chrom, pos) — ALT="." 表示 ref call
            if not matched:
                ref_key = (chrom, pos)
                if ref_key in ref_data:
                    ref_dp, seq_info, seq_fmt, seq_sample, seq_qual = ref_data[ref_key]
                    seq_matched_data = (seq_info, seq_fmt, seq_sample, seq_qual)
                    n_matched_ref += 1
                    matched = True
                    n_unchanged += 1

            if not matched:
                n_not_found += 1

            # 确保输出中所有基因型都是 0/0 或实际变异, 不出现缺失
            if fields[9] in ("./.", ".|.", "0|0"):
                fields[9] = "0/0"

            fout.write("\t".join(fields) + "\n")

            # 写 annotated 格式 TSV — 只输出发生变异的位点
            # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SM
            if gt_updated and seq_matched_data is not None:
                seq_info, seq_fmt, seq_sample, seq_qual = seq_matched_data
                # INFO: PX=GENE 前缀 + 测序VCF的完整INFO
                ann_info = f"PX={gene};{seq_info}" if gene != "." else seq_info
                # SM: 用测序数据的sample, 但GT替换为更新后的值
                sm_parts = seq_sample.split(":")
                fmt_fields_list = seq_fmt.split(":")
                if "GT" in fmt_fields_list:
                    gt_idx = fmt_fields_list.index("GT")
                    sm_parts[gt_idx] = fields[9]  # 使用更新后的GT
                ann_sm = ":".join(sm_parts)
                fann.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{seq_qual}\t{fields[6]}\t{ann_info}\t{seq_fmt}\t{ann_sm}\n")

    # ---- 报告 ----
    print("=" * 60)
    print("更新结果")
    print("=" * 60)
    print(f"  模板位点总数:        {n_total}")
    print(f"  匹配变异位点:        {n_matched_variant}")
    print(f"  匹配参考纯合位点:    {n_matched_ref}")
    print(f"  总匹配:              {n_matched_variant + n_matched_ref}")
    print(f"  基因型已更新:        {n_updated}  ({n_updated/n_total*100:.1f}%)")
    print(f"  基因型未变:          {n_unchanged}")
    print(f"  未匹配到测序数据:    {n_not_found}")

    if changed_records:
        # 按基因统计
        gene_changed = {}
        for rec in changed_records:
            g = rec[5]
            gene_changed[g] = gene_changed.get(g, 0) + 1
        print(f"\n  变更记录 (按基因):")
        for g, cnt in sorted(gene_changed.items(), key=lambda x: -x[1]):
            print(f"    {g}: {cnt}")

        # 按基因型统计
        gt_changed = {}
        for rec in changed_records:
            new_gt = rec[7]
            gt_changed[new_gt] = gt_changed.get(new_gt, 0) + 1
        print(f"\n  变更后的基因型分布:")
        for gt, cnt in sorted(gt_changed.items(), key=lambda x: -x[1]):
            print(f"    {gt}: {cnt}")

    print(f"\n  输出文件:")
    print(f"    更新后 VCF:    {OUT_VCF}")
    print(f"    变更记录:      {OUT_CHANGED}")
    print(f"    Annotated TSV: {OUT_ANNOTATED}")

def updata_outcall(MTRNR1_VAR):
    
    fout_fp = open(OUTCALL_TSV, "w", encoding="utf-8")
    
    if 'm.1095T>C' in MTRNR1_VAR:
        fout_fp.write('MT-RNR1\tm.1095T>C')
    elif 'm.1494C>T' in MTRNR1_VAR:
        fout_fp.write('MT-RNR1\tm.1494C>T')    
    elif 'm.1555A>G' in MTRNR1_VAR:
        fout_fp.write('MT-RNR1\tm.1555A>G')
    elif 'm.827A>G' in MTRNR1_VAR:
        fout_fp.write('MT-RNR1\tm.827A>G')
    else:
        fout_fp.write('MT-RNR1\tReference')
    return
    

# ======================== Step 3: NoCall 位点检测 ========================
def _read_bed(bed_path):
    """
    读取 BED 文件, 返回 {chrom: [(start, end), ...]}
    BED 坐标: 0-based, [start, end)
    """
    regions = {}
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((start, end))
    # 排序
    for chrom in regions:
        regions[chrom].sort()
    return regions


def _in_pcr_target(regions, chrom, pos):
    """
    检查 pos (1-based) 是否落在 PCR target 区域内。
    BED 是 0-based, VCF 是 1-based, 所以 pos-1 对应 BED 坐标。
    """
    if chrom not in regions:
        return False
    p = pos - 1  # 转成 0-based
    for start, end in regions[chrom]:
        if p >= start and p < end:
            return True
        if start > p:  # 已排序, 后面不会匹配
            break
    return False


def find_nocall_positions():
    """
    检测监测位点的 NoCall 状态。
    前提: 位点在 PharmCAT 模板中 且 在 PCR target 区域内。
    NoCall 条件 (满足任一):
      1. 位点在测序VCF中未出现 (未检测到)
      2. DP < DP_THRESHOLD (深度不够)
      3. QUAL < NOCALL_QUAL_MIN (质量太低)
    
    输出: {sample}.nocall.tsv
    """
    print("=" * 60)
    print("Step 3: NoCall 位点检测")
    print("=" * 60)
    
    # ---- 3a: 读取 PCR target BED ----
    pcr_regions = _read_bed(PCR_TARGET)
    n_pcr_regions = sum(len(v) for v in pcr_regions.values())
    print(f"  PCR target 区域数: {n_pcr_regions}")
    
    # ---- 3b: 扫描测序VCF, 记录每个位点的 DP 和 QUAL ----
    seq_positions = {}
    
    with open(SEQ_VCF, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            qual_str = fields[5]
            info = fields[7]
            fmt = fields[8]
            sample = fields[9] if len(fields) > 9 else ""
            
            # 跳过 impute 位点
            is_imputed = ("RAF=" in info) or ("DS" in fmt.split(":"))
            if is_imputed:
                continue
            
            # 解析 QUAL
            try:
                qual = float(qual_str) if qual_str != "." else 0.0
            except ValueError:
                qual = 0.0
            
            # 解析 DP
            fmt_fields = fmt.split(":")
            sample_fields = sample.split(":")
            fmt_idx = {ff: j for j, ff in enumerate(fmt_fields)}
            dp = 0
            if "DP" in fmt_idx:
                try:
                    dp = int(sample_fields[fmt_idx["DP"]])
                except (ValueError, IndexError):
                    dp = 0
            
            key = (chrom, pos)
            if key not in seq_positions or dp > seq_positions[key][0]:
                seq_positions[key] = (dp, qual)
    
    # ---- 3c: 读取模板, 逐个位点判断 NoCall ----
    nocall_records = []
    n_template = 0
    n_in_pcr = 0
    
    with open(TEMPLATE_VCF, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            rsid = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # 提取基因名
            gene = "."
            for kv in info.split(";"):
                if kv.startswith("PX="):
                    gene = kv[3:]
                    break
            
            n_template += 1
            
            # 检查是否在 PCR target 区域内
            if not _in_pcr_target(pcr_regions, chrom, pos):
                continue  # 不在 PCR target 内, 跳过
            n_in_pcr += 1
            
            key = (chrom, pos)
            
            # 判断 NoCall 原因
            reasons = []
            if key not in seq_positions:
                reasons.append("NOT_DETECTED")
            else:
                dp, qual = seq_positions[key]
                if dp < DP_THRESHOLD:
                    reasons.append(f"DP={dp}<{DP_THRESHOLD}")
                # QUAL 为 "." 表示参考纯合 (无变异), 不参与 QUAL 过滤
                # 只有实际 QUAL 值 < 阈值时才标记
                if qual > 0 and qual < NOCALL_QUAL_MIN:
                    reasons.append(f"QUAL={qual:.1f}<{NOCALL_QUAL_MIN:.0f}")
            
            if reasons:
                nocall_records.append((chrom, pos, rsid, ref, alt, gene, ";".join(reasons)))
    
    # ---- 3c: 写入 NoCall 文件 ----
    with open(OUT_NOCALL, "w", encoding="utf-8") as fout:
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tGENE\tREASON\n")
        for rec in nocall_records:
            fout.write("\t".join(str(x) for x in rec) + "\n")
    
    # ---- 报告 ----
    n_called = n_in_pcr - len(nocall_records)
    print(f"  模板位点总数:      {n_template}")
    print(f"  PCR target 内位点: {n_in_pcr}")
    print(f"  测序覆盖位点数:    {len(seq_positions)}")
    print(f"  NoCall 位点数:     {len(nocall_records)} ({len(nocall_records)/max(n_in_pcr,1)*100:.1f}%)")
    print(f"  有效 Call 位点:    {n_called} ({n_called/max(n_in_pcr,1)*100:.1f}%)")
    
    if nocall_records:
        # 按原因统计
        reason_count = {}
        for rec in nocall_records:
            for r in rec[6].split(";"):
                reason_count[r] = reason_count.get(r, 0) + 1
        print(f"\n  NoCall 原因统计:")
        for reason, cnt in sorted(reason_count.items(), key=lambda x: -x[1]):
            print(f"    {reason}: {cnt}")
        
        # 按基因统计
        gene_count = {}
        for rec in nocall_records:
            g = rec[5]
            gene_count[g] = gene_count.get(g, 0) + 1
        print(f"\n  NoCall 按基因统计:")
        for g, cnt in sorted(gene_count.items(), key=lambda x: -x[1]):
            print(f"    {g}: {cnt}")
    
    print(f"\n  NoCall 文件: {OUT_NOCALL}")
    print()
    return nocall_records


# ======================== Main ========================
def main():
    print(f"样本: {SAMPLE_NAME}")
    print(f"DP 阈值: >= {DP_THRESHOLD}")
    print(f"假阳性过滤: {'开启' if FP_FILTER_ENABLED else '关闭 (--no-fp-filter)'}")
    print()

    t0 = time.time()

    variant_data, ref_data = read_sequenced_vcf()
    update_template(variant_data, ref_data)
    updata_outcall(MTRNR1_VAR)
    find_nocall_positions()
    elapsed = time.time() - t0
    print(f"\n总耗时: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
