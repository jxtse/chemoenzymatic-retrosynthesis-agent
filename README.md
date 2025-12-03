# BRENDA SOAP Client

用法示例：
1. 安装依赖：`pip install -e .`
2. 在 `.env` 中设置账号：`BRENDA_EMAIL`（或 `BRENDA_EMIAL`）和 `BRENDA_PASSWORD`（明文密码，脚本内部会自动做 SHA-256）。
3. 调用示例：
   - 查询 Km：`python main.py --mode km --ec 1.1.1.1 --organism "*" --substrate "*"`
   - 查询反应式：`python main.py --mode reaction --ec 1.1.1.1 --organism "*" --reaction "*"`

脚本会自动读取 `.env`，通过 BRENDA SOAP 接口返回对应记录，并逐行输出。

## 统一的反应 JSON 导出

把本地的三个数据源转成一行一条的 LLM-ready JSONL：

- EnzyExtract：`python reaction_unifier.py --dataset enzyextract --input EnzyExtractDB_176463.parquet --output enzyextract.jsonl`
- BRENDA kcat/km：`python reaction_unifier.py --dataset brenda --input brenda_kcat_v3.parquet --output brenda.jsonl`
- BKMS 反应：`python reaction_unifier.py --dataset bkms --input Reactions_BKMS.csv --output bkms.jsonl`

可选：加 `--limit 1000` 快速抽样。

## BKMS 反应核心 + 动力学对齐

以 BKMS 反应为“主项”，按 EC Number 自动挂接 EnzyExtract/BRENDA 的动力学测量：

- `python reaction_merger.py --bkms Reactions_BKMS.csv --enzyextract EnzyExtractDB_176463.parquet --brenda brenda_kcat_v3.parquet --output merged.jsonl`

同样可以加 `--limit N` 先抽样。

输出视图：
- `--view full`（默认）：完整嵌套。
- `--view minimal`：保留方程、底物/产物名称、EC、测量的简要字段，附带 `substrate_ref/index`、标准化的 `kcat_std`（s^-1）、`km_std_mM` 和 log10。
- `--view flat`：一行一个测量样本，带 reaction 上下文，可直接做训练表。
