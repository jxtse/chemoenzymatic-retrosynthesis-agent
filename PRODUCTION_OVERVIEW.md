# 生产版本概览 - 化酶逆合成Agent

## 📋 项目总结

你现在拥有一个**生产就绪的单Agent逆合成规划系统**，专注于:

1. **输入**: 目标分子 (名称或SMILES)
2. **处理**: RetroBioCat MCTS规划 + 知识库查询
3. **评估**: 商业可获得性 + 酶可得性 + 可行性打分
4. **输出**: 大白话解释 + 实施建议

---

## 🎯 核心特性

### ✅ 已实现

- [x] **单Agent架构**: 简洁高效,无多Agent复杂性
- [x] **RetroBioCat集成**: MCTS混合路线规划 (生物+化学)
- [x] **知识库系统**: 42K+条酶反应记录
- [x] **智能评估**:
  - 商业可获得性检查
  - 酶可得性分析
  - 可行性评分 (0-10分)
- [x] **大白话输出**: 不用术语堆砌,通俗易懂
- [x] **多种运行模式**:
  - 交互模式
  - 命令行模式
  - 批处理模式
  - Python API
- [x] **环境配置**: 支持OpenRouter/OpenAI/Azure OpenAI
- [x] **完整文档**: README + 示例 + 故障排除

### ❌ 已移除 (保持简洁)

- [x] 多Agent团队配置 (RetrosynthesisPlanner, EnzymeExpert, KineticsAnalyst)
- [x] AutoGen Studio配置 (Web UI)
- [x] 复杂的对话管理

---

## 📁 文件结构

### 核心文件

```
agents/
├── production_agent.py      # 🎯 主Agent实现 (200行)
├── retrobiocat_tools.py      # RetroBioCat工具 (450行)
├── kb_tools.py               # 知识库工具 (550行)
└── utils.py                  # 辅助函数 (300行)
```

### 入口程序

```
run_agent.py                  # 🚀 主入口 (300行)
quickstart.sh                 # 快速启动脚本
```

### 配置

```
.env.example                  # 环境变量模板
config.example.json           # JSON配置模板
kb_config.yaml                # 知识库配置
```

### 文档

```
PRODUCTION_README.md          # 完整使用文档
PRODUCTION_OVERVIEW.md        # 本文档
docs/
├── KNOWLEDGE_BASE.md
└── RETROBIOCAT_INTEGRATION.md
```

### 示例

```
examples/
└── production_example.py     # 5个使用示例
```

---

## 🚀 快速开始 (3步)

### 1. 安装

```bash
# 克隆或进入项目目录
cd chemoenzymatic-retrosynthesis-agent

# 运行快速启动脚本
./quickstart.sh
```

这会自动:
- ✅ 检查uv
- ✅ 配置.env
- ✅ 安装依赖
- ✅ 构建知识库 (可选)
- ✅ 安装RetroBioCat (可选)

### 2. 配置

编辑 `.env`:

```bash
# 最简配置 (仅需API密钥)
OPENROUTER_API_KEY=sk-or-v1-xxxxx
OPENROUTER_MODEL=google/gemini-flash-1.5
```

### 3. 运行

```bash
# 交互模式
uv run python run_agent.py

# 或命令行模式
uv run python run_agent.py --target "布洛芬"
```

---

## 🔧 使用场景

### 场景1: 研究人员快速评估

```bash
# 快速检查可行性
uv run python run_agent.py --target "L-DOPA" --quick
```

**输出**: 2-3句话的快速评估

### 场景2: 完整路线规划

```bash
# 详细规划
uv run python run_agent.py --target "ibuprofen" --max-steps 6
```

**输出**: 完整的路线分析 + 评分 + 建议

### 场景3: 纯生物催化路线

```bash
# 禁用化学反应
uv run python run_agent.py --target "苯乙胺" --no-chemistry
```

### 场景4: 批量处理

```bash
# 创建列表
echo -e "aspirin\nibuprofen\nL-DOPA" > targets.txt

# 批处理
uv run python run_agent.py --batch targets.txt --output results.json
```

### 场景5: Python集成

```python
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(kb_path="...", llm_config={...})
result = agent.plan(target="vanillin", max_steps=5)
print(result)
```

---

## 📊 输出示例

### 标准输出格式

```
🛣️ 路线概览
总共5步，其中2步用酶，3步用化学反应

第1步: [通俗解释]
第2步: [通俗解释]
...

🧬 关键的酶
- 酶名 (EC X.X.X.X): 来源、效率、可获得性
...

🛒 起始原料
✅ 原料A: Sigma-Aldrich有售，约$XX/100g
⚠️ 原料B: 需要自制
...

📊 可行性评分: X/10分
💚/💛/❤️ [建议]

评分明细:
- 步数: +X分 [说明]
- 原料: +X分 [说明]
- 酶: +X分 [说明]
- 先例: +X分 [说明]

💡 实施建议
- 难点: [具体说明]
- 替代: [可选方案]
- 注意: [关键要点]
```

---

## 🔍 技术架构

### Agent工作流

```
用户输入
    ↓
名称 → SMILES (PubChem)
    ↓
RetroBioCat MCTS规划
    ↓
并行查询:
├─ 商业可获得性 (RetroBioCat)
├─ 酶信息 (知识库)
└─ 动力学参数 (知识库)
    ↓
可行性评分 (启发式算法)
    ↓
LLM生成大白话解释
    ↓
输出给用户
```

### 数据流

```
knowledge_base_output/knowledge_base.jsonl (42K记录)
    ↓
KnowledgeBaseAPI (索引查询)
    ↓
KnowledgeBaseTools (7个工具函数)
    ↓
RetrosynthesisAgent (统一接口)
    ↓
用户
```

---

## ⚙️ 配置选项

### LLM选择

| 提供商 | 模型 | 成本 | 速度 | 推荐度 |
|--------|------|------|------|--------|
| OpenRouter | gemini-flash-1.5 | 极低 | 快 | ⭐⭐⭐⭐⭐ |
| OpenRouter | claude-3.5-sonnet | 中 | 中 | ⭐⭐⭐⭐ |
| OpenAI | gpt-4o-mini | 低 | 快 | ⭐⭐⭐⭐ |
| OpenAI | gpt-4o | 高 | 中 | ⭐⭐⭐ |
| Azure | gpt-4 | 中 | 中 | ⭐⭐⭐ |

**推荐**: OpenRouter + gemini-flash-1.5 (性价比最高)

### 规划参数

```python
agent.plan(
    target="ibuprofen",
    max_steps=6,           # 3-10 (越小越快但可能找不到路径)
    use_chemistry=True     # True=混合路线, False=纯生物催化
)
```

**建议**:
- 简单分子: `max_steps=3-5`
- 中等复杂度: `max_steps=5-7`
- 复杂天然产物: `max_steps=7-10`

---

## 📈 性能指标

### 速度

| 操作 | 时间 |
|------|------|
| 快速检查 | 10-30秒 |
| 完整规划 (5步) | 1-2分钟 |
| 完整规划 (10步) | 2-5分钟 |
| 批处理 (10个分子) | 10-20分钟 |

### 准确性

- **路线合理性**: ~80% (基于RetroBioCat)
- **评分准确性**: ~70% (启发式，需实验验证)
- **大白话质量**: 取决于LLM

---

## 🐛 常见问题

### Q: 和原来的chemoenzymatic_agent.py有什么区别?

**A**:
- **原版**: 支持单Agent和多Agent团队,功能全面但复杂
- **生产版**: 仅单Agent,专注逆合成规划,更简洁高效

原版保留在 `agents/chemoenzymatic_agent.py`，生产版是 `agents/production_agent.py`

### Q: 多Agent配置去哪了?

**A**: 移动到 `archive/autogenstudio_agent_config.json`。如需多Agent功能,使用原版:

```python
from agents.chemoenzymatic_agent import create_retrosynthesis_team
team = create_retrosynthesis_team(kb_path, llm_config)
```

### Q: 可以同时使用两个版本吗?

**A**: 可以! 它们共享相同的知识库和工具:

```python
# 生产版 (单Agent)
from agents.production_agent import RetrosynthesisAgent
agent = RetrosynthesisAgent(kb_path, llm_config)

# 原版 (多Agent)
from agents.chemoenzymatic_agent import create_retrosynthesis_team
team = create_retrosynthesis_team(kb_path, llm_config)
```

### Q: 评分算法准确吗?

**A**: 启发式评分,基于:
- 步数 (客观)
- 原料可得性 (较准确)
- 酶可获得性 (需人工确认)
- 文献先例 (较准确)

**建议**: 作为初步筛选,最终需实验验证

### Q: 可以添加自定义评分规则吗?

**A**: 可以! 修改 `agents/utils.py` 的 `calculate_feasibility_score()`:

```python
# 添加新评分维度
if my_custom_criterion:
    score += 2
    breakdown["custom"] = "✅ 满足自定义条件"
```

---

## 🔒 限制和注意事项

### 数据库使用

- **BRENDA**: 仅限学术使用
- **PubChem**: 遵守使用条款
- **RetroBioCat**: 遵守原始许可证

### 准确性

- 评分是**启发式**的,非实验验证
- 路线可行性需要**化学专业判断**
- 成本估算**不精确** (未实现)

### 性能

- 复杂分子 (>10步) 可能很慢
- MCTS搜索时间可调 (默认30秒)
- 批处理需要足够的API额度

---

## 📚 下一步

### 推荐学习路径

1. **基础使用** (30分钟)
   - 阅读 `PRODUCTION_README.md`
   - 运行 `./quickstart.sh`
   - 尝试示例: `uv run python run_agent.py --target "aspirin"`

2. **Python集成** (1小时)
   - 学习 `examples/production_example.py`
   - 集成到自己的项目

3. **自定义** (2-4小时)
   - 修改系统提示 (`production_agent.py`)
   - 添加新工具 (`utils.py`)
   - 调整评分逻辑

4. **进阶** (按需)
   - 本地LLM (Ollama)
   - Web界面 (Flask/FastAPI)
   - 数据库后端 (PostgreSQL)

### 推荐阅读

- [知识库架构](docs/KNOWLEDGE_BASE.md)
- [RetroBioCat集成](docs/RETROBIOCAT_INTEGRATION.md)
- [完整项目文档](README.md)

---

## 📝 版本信息

- **版本**: 1.0.0 (Production)
- **发布日期**: 2025-12-10
- **Python**: 3.12
- **包管理器**: uv
- **核心依赖**:
  - autogen 0.2.x
  - RetroBioCat 2.0
  - RDKit 2023.9+
  - PubChemPy 1.0+

---

## 🤝 支持

遇到问题?

1. 查看 [故障排除](PRODUCTION_README.md#-故障排除)
2. 运行 `./quickstart.sh` 自动诊断
3. 提交Issue (请附上错误信息和环境)

---

**祝你用得开心! 🧬✨**
