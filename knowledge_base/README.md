# 统一化学-酶-动力学知识库

一个整合 **BKMS, BRENDA, EnzyExtract, KEGG, USPTO, UniProt, PubChem, ZINC, RetroBioCat** 等多个数据库的统一知识层,提供标准化的 JSON API 接口。

## 核心特性

✅ **9个数据源整合** - 本地文件(BKMS, BRENDA, EnzyExtract, RetroBioCat) + REST API(KEGG, UniProt, PubChem) + SOAP API(BRENDA)
✅ **统一Schema** - 所有数据转换为标准JSON格式
✅ **智能缓存** - 内存+磁盘双层缓存,减少API调用
✅ **速率控制** - 自动遵守各API速率限制
✅ **灵活查询** - 按EC号、化合物、反应、动力学参数查询
✅ **多种输出** - JSONL, Parquet
✅ **CLI + Python API** - 命令行工具 + 编程接口

## 快速开始

### 安装

```bash
pip install pandas pyarrow requests zeep pyyaml tqdm
```

### 5分钟上手

```bash
# 1. 创建配置
python -m knowledge_base.cli init-config

# 2. 编辑 kb_config.yaml 设置数据路径

# 3. 构建知识库
python -m knowledge_base.cli build --config kb_config.yaml

# 4. 查询
python -m knowledge_base.cli query \
    --kb knowledge_base_output/knowledge_base.jsonl \
    --ec 1.1.1.1 \
    --output result.json
```

### Python API

```python
from knowledge_base import KnowledgeBaseAPI

# 加载
api = KnowledgeBaseAPI("knowledge_base_output/knowledge_base.jsonl")
api.load()

# 查询
result = api.query_by_ec("1.1.1.1")
print(f"Found {result['count']} records")

# 动力学参数
kinetics = api.get_kinetics(ec_number="1.1.1.1", organism="E. coli")
for k in kinetics['kinetics']:
    print(f"kcat={k['kcat']} {k['kcat_unit']}, Km={k['km']} {k['km_unit']}")
```

## 架构

```
数据源 → Connectors → Unified Schema → Builder → Storage → API → JSON输出
```

- **Connectors**: 9个数据库适配器,统一查询接口
- **Unified Schema**: 标准化数据结构(EC号、酶、反应、动力学、条件)
- **Builder**: 多源数据整合、索引构建、并行处理
- **API**: EC查询、化合物查询、反应查询、动力学提取

## 主要组件

| 组件 | 文件 | 功能 |
|------|------|------|
| **配置** | `config.py` | 数据源配置、路径管理 |
| **连接器** | `connectors/*.py` | 9个数据库适配器 |
| **Schema** | `schema.py` | 统一数据结构定义 |
| **构建器** | `builder.py` | 知识库构建流程 |
| **API** | `api.py` | 统一查询接口 |
| **CLI** | `cli.py` | 命令行工具 |

## 数据源支持

| 数据库 | 类型 | 内容 | 状态 |
|--------|------|------|------|
| BKMS | 本地 | 反应方程、EC号、通路 | ✅ 完整 |
| BRENDA | 本地+API | kcat, Km, 生物体 | ✅ 完整 |
| EnzyExtract | 本地 | 动力学、序列、SMILES | ✅ 完整 |
| KEGG | API | 通路、酶、化合物 | ✅ 完整 |
| UniProt | API | 蛋白序列、结构 | ✅ 完整 |
| PubChem | API | 化合物性质、结构 | ✅ 完整 |
| RetroBioCat | 本地 | 反应模板、酶类型 | ✅ 完整 |
| ZINC | API | 可购买化合物 | ⚠️ 存根 |
| USPTO | API | 专利数据 | ⚠️ 存根 |

## 查询示例

### CLI

```bash
# EC号查询
python -m knowledge_base.cli query --kb kb.jsonl --ec 1.1.1.1

# 化合物查询
python -m knowledge_base.cli query --kb kb.jsonl --compound glucose --role substrate

# 反应查询
python -m knowledge_base.cli query --kb kb.jsonl --reaction \
    --substrates "ATP,glucose" --products "ADP,glucose-6-phosphate"

# 动力学参数
python -m knowledge_base.cli query --kb kb.jsonl --kinetics \
    --ec-filter 1.1.1.1 --organism-filter "Escherichia coli"

# 统计
python -m knowledge_base.cli query --kb kb.jsonl --stats
```

### Python

```python
from knowledge_base import KnowledgeBaseAPI

api = KnowledgeBaseAPI("kb.jsonl")
api.load()

# EC查询
result = api.query_by_ec("1.1.1.1", include_sources=["BRENDA", "BKMS"])

# 化合物查询
result = api.query_by_compound("glucose", role="substrate")

# 反应查询
result = api.query_by_reaction(
    substrates=["ATP", "glucose"],
    products=["ADP", "glucose-6-phosphate"]
)

# 动力学查询
kinetics = api.get_kinetics(
    ec_number="1.1.1.1",
    substrate="ethanol",
    organism="Saccharomyces cerevisiae"
)

# 统计
stats = api.get_statistics()
```

## 配置示例

```yaml
# kb_config.yaml
data_dir: "."
output_dir: "knowledge_base_output"
output_format: "jsonl"  # or "parquet"
compress_output: false

bkms:
  enabled: true
  path: "Reactions_BKMS.csv"

brenda:
  enabled: true
  path: "brenda_kcat_v3.parquet"
  api_url: "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
  rate_limit: 0.5

kegg:
  enabled: true
  rate_limit: 0.2  # 5 sec between requests

uniprot:
  enabled: true
  rate_limit: 1.0

pubchem:
  enabled: true
  rate_limit: 5.0
```

## 统一 Schema

```json
{
  "id": "BKMS:0",
  "source": {"dataset": "BKMS", "row_index": 0},
  "ec_numbers": ["1.1.1.1"],
  "primary_ec": "1.1.1.1",
  "enzyme": {
    "name": "alcohol dehydrogenase",
    "uniprot_ids": ["P00330"],
    "sequence": "MSKL...",
    "organism": "Escherichia coli"
  },
  "reaction": {
    "equation_text": "ethanol + NAD+ => acetaldehyde + NADH",
    "substrates": [{"name": "ethanol", "smiles": "CCO", "cid": "702"}],
    "products": [{"name": "acetaldehyde"}]
  },
  "kinetics": {
    "kcat": {"value": 150.0, "unit": "s^-1"},
    "km": {"value": 0.5, "unit": "mM"}
  },
  "conditions": {
    "pH": 7.4,
    "temperature": 37.0,
    "buffer": "phosphate"
  }
}
```

## 性能

- **本地数据加载**: 2-5秒 (230K+ 记录)
- **索引构建**: ~1秒
- **查询响应**: <10ms (本地)
- **API丰富化**: 0.2-5秒/EC号

## 环境变量

```bash
export BRENDA_EMAIL="your.email@example.com"
export BRENDA_PASSWORD="your_password"
```

## 文档

详细文档: [KNOWLEDGE_BASE.md](../KNOWLEDGE_BASE.md)

## 示例

```bash
python examples/build_and_query_kb.py
```

## License

遵循各数据源许可证
