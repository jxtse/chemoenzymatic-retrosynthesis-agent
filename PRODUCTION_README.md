# 化酶逆合成规划Agent - 生产版本

生产就绪的单Agent实现，专注于**逆合成路线规划 + 可行性评估**。

## 🎯 核心功能

输入目标分子 → RetroBioCat规划路线 → 评估可行性 → 大白话输出建议

### 工作流程

1. **理解任务**: 分子名称 → SMILES转换
2. **规划路线**: RetroBioCat MCTS混合规划 (生物酶 + 化学)
3. **评估路径**:
   - ✅ 检查商业可获得性
   - 🧬 查询酶的来源和效率
   - 📊 计算可行性评分 (0-10分)
4. **大白话输出**: 不用术语堆砌，像给高中生讲课一样清楚

---

## 🚀 快速开始

### 1. 安装依赖

```bash
# 确保使用uv作为包管理器
uv sync

# 安装RetroBioCat (必需)
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### 2. 配置环境变量

创建 `.env` 文件:

```bash
cp .env.example .env
```

编辑 `.env` 并设置API密钥:

```env
# 推荐: OpenRouter (支持多种模型)
OPENROUTER_API_KEY=sk-or-v1-xxxxx
OPENROUTER_MODEL=google/gemini-flash-1.5

# 或使用 OpenAI
# OPENAI_API_KEY=sk-xxxxx
# OPENAI_MODEL=gpt-4o-mini

# 或使用 Azure OpenAI
# AZURE_OPENAI_API_KEY=xxxxx
# AZURE_OPENAI_ENDPOINT=https://xxx.openai.azure.com/
# AZURE_OPENAI_DEPLOYMENT=gpt-4
```

### 3. 构建知识库 (首次运行)

```bash
uv run python -m knowledge_base.cli build --config kb_config.yaml
```

这会创建 `knowledge_base_output/knowledge_base.jsonl` (~76MB, 42K+条记录)

### 4. 运行Agent

#### 交互模式 (推荐)

```bash
uv run python run_agent.py
```

然后输入目标分子:
```
请输入目标分子: 布洛芬
```

#### 命令行模式

```bash
# 完整规划
uv run python run_agent.py --target "ibuprofen"

# 快速检查
uv run python run_agent.py --target "L-DOPA" --quick

# 纯生物催化路线
uv run python run_agent.py --target "苯乙胺" --no-chemistry

# 指定最大步数
uv run python run_agent.py --target "aspirin" --max-steps 5
```

#### 批处理模式

```bash
# 创建目标列表文件
echo -e "ibuprofen\nL-DOPA\naspirin" > targets.txt

# 批量处理
uv run python run_agent.py --batch targets.txt --output results.json
```

---

## 📁 项目结构 (生产版本)

```
chemoenzymatic-retrosynthesis-agent/
├── agents/
│   ├── production_agent.py       # 🎯 生产Agent (主要实现)
│   ├── retrobiocat_tools.py       # RetroBioCat工具
│   ├── kb_tools.py                # 知识库工具
│   └── utils.py                   # 辅助函数 (name→SMILES, 评分等)
│
├── knowledge_base/                # 知识库系统
│   ├── api.py                     # 查询API
│   ├── builder.py                 # 构建器
│   └── connectors/                # 9个数据库连接器
│
├── knowledge_base_output/         # 构建好的知识库
│   └── knowledge_base.jsonl       # 42K+条记录
│
├── run_agent.py                   # 🚀 主入口程序
├── config.example.json            # 配置示例
├── .env.example                   # 环境变量示例
│
├── examples/
│   └── production_example.py      # 使用示例
│
├── docs/                          # 文档
│   ├── KNOWLEDGE_BASE.md
│   ├── RETROBIOCAT_INTEGRATION.md
│   └── PRODUCTION_README.md       # 本文档
│
└── datasets/                      # 数据源
    ├── Reactions_BKMS.csv         # 11MB
    ├── brenda_kcat_v3.parquet     # 3.4MB
    └── EnzyExtractDB_176463.parquet # 9.7MB
```

---

## 🔧 核心组件

### 1. RetrosynthesisAgent

主Agent类，提供两个核心方法:

```python
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config=llm_config
)

# 完整规划
result = agent.plan(
    target="ibuprofen",
    max_steps=6,
    use_chemistry=True
)

# 快速检查
quick_result = agent.quick_check("L-DOPA")
```

### 2. 工具函数

**分子处理**:
- `name_to_smiles(compound_name)`: 名称 → SMILES (PubChem)
- `analyze_molecule_properties(smiles)`: 分析分子性质

**路线规划**:
- `plan_biocatalytic_route()`: MCTS混合路线规划
- `find_enzymatic_reactions()`: 单步酶促反应
- `check_commercial_availability()`: 商业可获得性

**知识库查询**:
- `search_enzyme_by_ec()`: EC号查酶
- `get_kinetic_parameters()`: 查kcat/Km

**评估**:
- `calculate_feasibility_score()`: 可行性打分 (0-10)

---

## 📊 输出格式

Agent会用**大白话**输出，包含:

### 🛣️ 路线概览
- 总共X步，哪些用酶哪些用化学
- 每步通俗解释在干什么

### 🧬 关键的酶
- 需要哪些酶、从哪来 (大肠杆菌、酵母...)
- 效率怎么样 (kcat, Km)
- 好不好表达/获取

### 🛒 起始原料
- 哪些能买到
- 哪些需要自制
- 大概价格/可获得性

### 📊 可行性评分 (X/10分)

评分标准:
- **步数** (3步内+3分, 4-6步+2分, 7步以上+1分)
- **原料可得性** (全能买+3分, 部分+2分, 都不能+1分)
- **酶可获得性** (常见酶+3分, 需工程+2分, 罕见+1分)
- **文献先例** (有先例+2分)

解读:
- **8-10分**: 💚 强烈推荐！
- **5-7分**: 💛 可以试，注意难点
- **1-4分**: ❤️ 不推荐，换思路

### 💡 实施建议
- 难点在哪
- 替代方案
- 注意事项

---

## 🎨 使用示例

### 示例1: 布洛芬合成

```bash
uv run python run_agent.py --target "ibuprofen"
```

**预期输出** (示意):

```
🛣️ 路线概览
总共5步，其中2步用酶，3步用化学反应

第1步: 把异丁苯变成异丁苯酮 (化学氧化)
第2步: 酮变成醇 (用醇脱氢酶，EC 1.1.1.1)
第3步: ...

🧬 关键的酶
- 醇脱氢酶 (EC 1.1.1.1): 来自大肠杆菌，kcat=150/s，效率高
- 容易获得，Sigma有卖重组酶

🛒 起始原料
✅ 异丁苯: Sigma-Aldrich有售，约$50/100g
✅ NAD+辅因子: 商业可得

📊 可行性评分: 8/10分
💚 强烈推荐！这条路线很靠谱

评分明细:
- 步数合理 (5步) +2分
- 原料全能买到 +3分
- 酶常见易得 +3分
- 有文献先例 +0分

💡 实施建议
难点: 第2步需要辅因子循环系统
建议: 用葡萄糖脱氢酶再生NAD+
注意: 反应需要在pH 7.4, 30°C进行
```

### 示例2: 快速检查

```python
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(kb_path="knowledge_base_output/knowledge_base.jsonl")
result = agent.quick_check("L-DOPA")

# 输出: 简要可行性分析 (2-3句话)
```

### 示例3: Python API

```python
import os
from agents.production_agent import RetrosynthesisAgent

# 配置
llm_config = {
    "config_list": [{
        "model": "google/gemini-flash-1.5",
        "api_key": os.getenv("OPENROUTER_API_KEY"),
        "base_url": "https://openrouter.ai/api/v1",
    }],
    "temperature": 0.7,
}

# 创建Agent
agent = RetrosynthesisAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config=llm_config
)

# 规划
result = agent.plan(
    target="vanillin",  # 香草醛
    max_steps=5,
    use_chemistry=True
)

print(result)
```

---

## ⚙️ 配置选项

### LLM配置

支持多种LLM提供商:

**OpenRouter** (推荐):
```python
{
    "config_list": [{
        "model": "google/gemini-flash-1.5",  # 便宜快速
        "api_key": "sk-or-v1-xxxxx",
        "base_url": "https://openrouter.ai/api/v1",
    }],
    "temperature": 0.7,
}
```

**OpenAI**:
```python
{
    "config_list": [{
        "model": "gpt-4o-mini",
        "api_key": "sk-xxxxx",
    }],
}
```

**Azure OpenAI**:
```python
{
    "config_list": [{
        "model": "gpt-4",
        "api_type": "azure",
        "api_key": "xxxxx",
        "base_url": "https://xxx.openai.azure.com/",
        "api_version": "2024-02-01",
    }],
}
```

### RetroBioCat配置

```python
agent.plan(
    target="ibuprofen",
    max_steps=6,           # 最大步数 (3-10)
    use_chemistry=True     # True=混合路线, False=纯生物催化
)
```

---

## 🔍 知识库

### 数据来源

| 数据库 | 类型 | 记录数 | 状态 |
|--------|------|--------|------|
| BKMS | 本地CSV | 42,539 | ✅ |
| BRENDA | 本地Parquet | ~162 | ✅ |
| EnzyExtract | 本地Parquet | ~170K | ✅ |
| KEGG | REST API | - | ⚠️ 可选 |
| UniProt | REST API | - | ⚠️ 可选 |
| PubChem | REST API | - | ⚠️ 可选 |

**当前知识库**: 42,701条统一记录 (76MB JSONL)

### 重新构建知识库

```bash
# 编辑配置
vim kb_config.yaml

# 构建
uv run python -m knowledge_base.cli build --config kb_config.yaml

# 查询统计
uv run python -m knowledge_base.cli query --kb knowledge_base_output/knowledge_base.jsonl --stats
```

---

## 🐛 故障排除

### 问题1: RetroBioCat未安装

**症状**: `ImportError: No module named 'rbc2'`

**解决**:
```bash
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### 问题2: 知识库未找到

**症状**: `知识库文件不存在`

**解决**:
```bash
uv run python -m knowledge_base.cli build --config kb_config.yaml
```

### 问题3: API密钥未设置

**症状**: `未找到API密钥`

**解决**:
```bash
# 检查.env文件
cat .env

# 确保设置了以下之一:
# OPENROUTER_API_KEY
# OPENAI_API_KEY
# AZURE_OPENAI_API_KEY
```

### 问题4: RDKit导入错误

**症状**: `cannot import name 'Chem' from 'rdkit'`

**解决**:
```bash
uv add rdkit
```

---

## 📚 进阶功能

### 自定义系统提示

修改 `agents/production_agent.py` 的 `_get_system_message()` 方法来定制Agent行为。

### 添加新工具

在 `agents/utils.py` 添加函数，然后在 `production_agent.py` 的 `_register_tools()` 中注册:

```python
def my_custom_tool(arg: str) -> str:
    """自定义工具"""
    return json.dumps({"result": "..."})

# 在_register_tools中:
tools["my_custom_tool"] = my_custom_tool
```

### 集成到其他系统

```python
# 作为Python库使用
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(kb_path="...", llm_config={...})

# API模式
result = agent.plan(target="aspirin")

# 解析结果
if "💚" in result:
    print("高可行性路线!")
```

---

## 📖 相关文档

- [知识库架构](docs/KNOWLEDGE_BASE.md)
- [RetroBioCat集成](docs/RETROBIOCAT_INTEGRATION.md)
- [完整项目文档](README.md)

---

## 🙋 常见问题

**Q: 和原来的多Agent版本有什么区别?**

A: 生产版本是**单Agent**，专注于逆合成规划，流程更简洁高效。多Agent版本功能更全但复杂度更高。

**Q: 可以离线运行吗?**

A: 知识库可以离线，但LLM需要API调用。可以配置本地LLM (如Ollama)。

**Q: 支持哪些分子?**

A: 理论上支持所有有机小分子 (<1000 Da)。复杂天然产物可能需要更长时间规划。

**Q: 评分准确吗?**

A: 评分是启发式的，基于步数、原料、酶可获得性和文献先例。最终需要实验验证。

**Q: 可以商用吗?**

A: 代码本身可以，但需要注意:
- LLM API使用条款
- RetroBioCat许可证
- 数据库使用限制 (BRENDA学术使用only)

---

## 📝 许可证

本项目代码遵循 MIT License。

**注意**:
- BRENDA数据库仅限学术使用
- RetroBioCat遵循其原始许可证
- LLM API使用需遵循各提供商条款

---

## 🤝 贡献

欢迎提Issue和PR!

重点改进方向:
- 更准确的可行性评分
- 更多数据库集成
- 成本估算功能
- Web界面

---

**最后更新**: 2025-12-10
**版本**: 1.0.0 (Production)
