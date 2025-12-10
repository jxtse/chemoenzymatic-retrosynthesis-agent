# 快速参考卡 - 化酶逆合成Agent

## 🚀 5分钟快速开始

```bash
# 1. 快速启动
./quickstart.sh

# 2. 运行Agent
uv run python run_agent.py

# 3. 输入目标分子
请输入目标分子: 布洛芬
```

---

## 📝 常用命令

### 交互模式
```bash
uv run python run_agent.py
```

### 单分子规划
```bash
uv run python run_agent.py --target "ibuprofen"
```

### 快速检查
```bash
uv run python run_agent.py --target "L-DOPA" --quick
```

### 批处理
```bash
echo -e "aspirin\nibuprofen\nL-DOPA" > targets.txt
uv run python run_agent.py --batch targets.txt --output results.json
```

### 自定义参数
```bash
# 纯生物催化
uv run python run_agent.py --target "苯乙胺" --no-chemistry

# 限制步数
uv run python run_agent.py --target "vanillin" --max-steps 5
```

---

## 🔧 配置速查

### 最简配置 (.env)
```env
OPENROUTER_API_KEY=sk-or-v1-xxxxx
OPENROUTER_MODEL=google/gemini-flash-1.5
```

### Python API
```python
from agents.production_agent import RetrosynthesisAgent

agent = RetrosynthesisAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config={
        "config_list": [{
            "model": "google/gemini-flash-1.5",
            "api_key": "sk-or-v1-xxxxx",
            "base_url": "https://openrouter.ai/api/v1",
        }],
        "temperature": 0.7,
    }
)

# 完整规划
result = agent.plan(target="ibuprofen", max_steps=6)

# 快速检查
quick = agent.quick_check("L-DOPA")
```

---

## 📊 评分解读

| 分数 | 含义 | 建议 |
|------|------|------|
| 8-10 | 💚 高可行性 | 强烈推荐尝试 |
| 5-7  | 💛 中等可行性 | 可以试,注意难点 |
| 1-4  | ❤️ 低可行性 | 不推荐,换思路 |

**评分维度**:
- 步数 (0-3分)
- 原料可获得性 (0-3分)
- 酶可获得性 (0-2分)
- 文献先例 (0-2分)

---

## 🐛 故障速查

### 问题: `No module named 'rbc2'`
```bash
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git
```

### 问题: 知识库未找到
```bash
uv run python -m knowledge_base.cli build --config kb_config.yaml
```

### 问题: API密钥未设置
```bash
# 检查.env
cat .env

# 确保有以下之一
# OPENROUTER_API_KEY=...
# OPENAI_API_KEY=...
# AZURE_OPENAI_API_KEY=...
```

### 问题: RDKit错误
```bash
uv add rdkit
```

---

## 📁 文件速查

| 文件 | 用途 |
|------|------|
| `run_agent.py` | 主入口程序 |
| `agents/production_agent.py` | 生产Agent实现 |
| `agents/utils.py` | 工具函数 |
| `.env` | 环境变量配置 |
| `PRODUCTION_README.md` | 完整文档 |
| `quickstart.sh` | 快速启动脚本 |

---

## 🔗 相关链接

- **完整文档**: `PRODUCTION_README.md`
- **项目概览**: `PRODUCTION_OVERVIEW.md`
- **改动总结**: `PRODUCTION_CHANGES.md`
- **知识库**: `docs/KNOWLEDGE_BASE.md`
- **RetroBioCat**: `docs/RETROBIOCAT_INTEGRATION.md`

---

## 💡 最佳实践

### 1. 首次使用
```bash
./quickstart.sh  # 自动配置
```

### 2. 日常使用
```bash
uv run python run_agent.py  # 交互模式
```

### 3. 批量评估
```bash
# 创建列表 → 批处理 → 查看JSON结果
```

### 4. Python集成
```python
# 导入 → 配置 → plan() → 解析结果
```

---

## ⚙️ 参数速查

### plan()
```python
agent.plan(
    target="ibuprofen",      # 分子名称或SMILES
    max_steps=6,             # 3-10步
    use_chemistry=True       # True=混合, False=纯生物
)
```

### 命令行
```bash
--target "分子"        # 目标分子
--quick               # 快速检查
--batch file.txt      # 批处理
--output results.json # 输出文件
--max-steps 5         # 最大步数
--no-chemistry        # 纯生物催化
```

---

## 🎯 推荐工作流

### 研发筛选
```
目标列表 → 批处理 → 筛选高分(8+) → 详细规划
```

### 快速评估
```
分子名 → quick_check() → 2-3句话结论
```

### 深入分析
```
分子 → plan() → 完整报告 → 人工审查
```

---

保存此文档以便快速查阅！
