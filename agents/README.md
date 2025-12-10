# Autogen AI Agents for Chemoenzymatic Retrosynthesis

基于 Microsoft Autogen 的多智能体系统,用于化学酶促逆合成规划,集成统一知识库做检索增强生成 (RAG)。

## 概览

```
┌─────────────────────────────────────────────────────────┐
│           Unified Knowledge Base (230K+ records)         │
│   BKMS | BRENDA | EnzyExtract | KEGG | UniProt | ...   │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│              KnowledgeBaseTools (RAG Layer)              │
│  - search_enzyme_by_ec()                                │
│  - search_reactions_by_compound()                       │
│  - find_retrosynthesis_pathway()                        │
│  - get_kinetic_parameters()                             │
│  - get_enzyme_sequence()                                │
│  - compare_enzymes()                                    │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│                    Autogen Agents                        │
│                                                          │
│  Single Agent:                                           │
│  └─ ChemoenzymaticAgent (通用逆合成专家)                  │
│                                                          │
│  Multi-Agent Team:                                       │
│  ├─ RetrosynthesisPlanner (路径设计)                     │
│  ├─ EnzymeExpert (酶选择)                                │
│  └─ KineticsAnalyst (动力学分析)                         │
└─────────────────────────────────────────────────────────┘
                          ↓
                   User Applications
```

## 核心组件

### 1. KnowledgeBaseTools (`kb_tools.py`)

封装知识库 API 为 Autogen 工具函数,提供 7 个核心工具:

| 工具 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `search_enzyme_by_ec` | EC 号查询酶 | EC 号 | 酶信息、反应、动力学 |
| `search_reactions_by_compound` | 化合物查询反应 | 化合物名、角色 | 相关反应列表 |
| `find_retrosynthesis_pathway` | 查找逆合成路径 | 目标化合物、底物 | 可能的酶促路径 |
| `get_kinetic_parameters` | 获取动力学参数 | EC/底物/生物体 | kcat, Km 数据 |
| `get_enzyme_sequence` | 获取酶序列 | EC 号、生物体 | 氨基酸序列 |
| `compare_enzymes` | 比较多个酶 | EC 号列表 | 动力学对比 |
| `get_statistics` | 知识库统计 | - | 总体统计数据 |

**特性**:
- ✅ 自动格式化为 JSON (LLM 友好)
- ✅ OpenAI Function Calling 兼容
- ✅ 限制输出长度 (避免 token 溢出)
- ✅ 错误处理和日志

### 2. ChemoenzymaticAgent (`chemoenzymatic_agent.py`)

单一智能体,专注逆合成规划:

```python
from agents import ChemoenzymaticAgent

agent = ChemoenzymaticAgent(
    kb_path="knowledge_base.jsonl",
    llm_config={"model": "gpt-4", "api_key": "..."}
)

response = agent.chat(
    "Design an enzymatic route for glucose-6-phosphate from glucose and ATP"
)
```

**系统提示**:
- 化学酶促逆合成专家
- 使用知识库工具查询数据
- 提供 EC 号、动力学参数、生物体来源
- 考虑辅因子需求和生物体兼容性

### 3. Multi-Agent Team

多智能体协作系统:

```python
from agents import create_retrosynthesis_team

team = create_retrosynthesis_team(
    kb_path="knowledge_base.jsonl",
    llm_config={"model": "gpt-4", "api_key": "..."}
)

# 团队成员
# - RetrosynthesisPlanner: 设计合成路径
# - EnzymeExpert: 选择最优酶
# - KineticsAnalyst: 分析动力学
# - GroupChatManager: 协调对话

team["user_proxy"].initiate_chat(
    team["manager"],
    message="Design a pathway for L-DOPA synthesis"
)
```

**智能体分工**:

| 智能体 | 专长 | 职责 |
|--------|------|------|
| RetrosynthesisPlanner | 逆合成设计 | 分析目标、识别前体、设计多步路径 |
| EnzymeExpert | 酶选择 | 比较酶变体、评估底物特异性、推荐候选 |
| KineticsAnalyst | 反应工程 | 分析 kcat/Km、预测反应速率、识别限速步骤 |

### 4. RetrosynthesisSession

高层封装,简化使用:

```python
from agents import RetrosynthesisSession

session = RetrosynthesisSession(
    kb_path="knowledge_base.jsonl",
    llm_config={...},
    use_multi_agent=True  # 使用多智能体团队
)

# 规划合成
session.plan_synthesis(
    target_compound="L-DOPA",
    available_substrates=["tyrosine", "phenylalanine"],
    constraints={"organism": "E. coli"}
)

# 选择酶
session.select_enzyme(
    ec_number="1.14.16.2",
    criteria={"kcat_min": 10, "organism_preference": "E. coli"}
)

# 优化条件
session.optimize_conditions(
    ec_number="1.14.16.2",
    substrate="tyrosine"
)
```

## 安装

```bash
# 1. 安装 Autogen
pip install pyautogen

# 或安装完整依赖
pip install pyautogen openai python-dotenv

# 2. 设置 API 密钥
export OPENAI_API_KEY="your-openai-api-key"

# 或使用 .env 文件
echo "OPENAI_API_KEY=your-key" > .env
```

## 快速开始

### 1. 单智能体查询

```python
from agents import ChemoenzymaticAgent

# 创建智能体
agent = ChemoenzymaticAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config={
        "config_list": [{
            "model": "gpt-4",
            "api_key": "your-api-key",
        }],
        "temperature": 0.7,
    }
)

# 提问
response = agent.chat("""
Design an enzymatic synthesis route for glucose-6-phosphate.

Available: glucose, ATP
Organism: E. coli

Provide:
1. Reaction steps
2. EC numbers
3. Kinetic parameters
4. Alternative routes
""")
```

### 2. 多智能体协作

```python
from agents import create_retrosynthesis_team

# 创建团队
team = create_retrosynthesis_team(
    kb_path="knowledge_base.jsonl",
    llm_config={...}
)

# 复杂任务
team["user_proxy"].initiate_chat(
    team["manager"],
    message="""
Synthesize L-DOPA from tyrosine.

Requirements:
1. Multi-step pathway design (Planner)
2. Enzyme selection for each step (EnzymeExpert)
3. Kinetic analysis (KineticsAnalyst)
4. Overall efficiency evaluation
"""
)
```

### 3. 高层会话

```python
from agents import RetrosynthesisSession

session = RetrosynthesisSession(
    kb_path="knowledge_base.jsonl",
    llm_config={...},
    use_multi_agent=True
)

# 一行规划
session.plan_synthesis(
    target_compound="acetaldehyde",
    available_substrates=["ethanol", "NAD+"]
)
```

### 4. 直接使用工具 (无 LLM)

```python
from agents import KnowledgeBaseTools

tools = KnowledgeBaseTools("knowledge_base.jsonl")

# 查询酶
result = tools.search_enzyme_by_ec("1.1.1.1")
print(result)

# 查找路径
result = tools.find_retrosynthesis_pathway(
    target_compound="glucose-6-phosphate",
    available_substrates=["glucose", "ATP"]
)
print(result)

# 获取动力学
result = tools.get_kinetic_parameters(
    ec_number="1.1.1.1",
    organism="E. coli"
)
print(result)
```

## 示例

运行完整示例:

```bash
# 设置 API 密钥
export OPENAI_API_KEY="your-key"

# 运行示例
python examples/autogen_retrosynthesis.py
```

示例包括:
1. 单智能体逆合成规划
2. 多智能体团队协作
3. 高层会话接口
4. 直接工具使用 (无需 LLM)
5. 自定义智能体 (蛋白工程专家)

## 高级用法

### 自定义智能体

```python
from autogen import AssistantAgent
from agents import KnowledgeBaseTools

kb_tools = KnowledgeBaseTools("knowledge_base.jsonl")

# 创建专门的蛋白工程智能体
protein_engineer = AssistantAgent(
    name="ProteinEngineer",
    system_message="""You are a protein engineering expert.
Analyze enzyme sequences and suggest mutations for improved activity.""",
    llm_config={...}
)

# 注册工具
protein_engineer.register_function(
    function_map={
        "get_enzyme_sequence": kb_tools.get_enzyme_sequence,
        "compare_enzymes": kb_tools.compare_enzymes,
    }
)
```

### 集成外部工具

```python
from agents import create_retrosynthesis_team

team = create_retrosynthesis_team(...)

# 添加 RDKit 工具用于化学计算
def calculate_molecular_weight(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        return f"Molecular weight: {mw:.2f} g/mol"
    return "Invalid SMILES"

# 注册到智能体
team["planner"].register_function(
    function_map={"calculate_molecular_weight": calculate_molecular_weight}
)
```

### 保存对话历史

```python
from agents import ChemoenzymaticAgent
import json

agent = ChemoenzymaticAgent(...)

# 对话
agent.chat("Design route for compound X")

# 保存历史
chat_history = agent.agent.chat_messages
with open("conversation.json", "w") as f:
    json.dump(chat_history, f, indent=2)
```

## LLM 配置

### OpenAI

```python
llm_config = {
    "config_list": [{
        "model": "gpt-4",  # 或 gpt-3.5-turbo
        "api_key": os.environ["OPENAI_API_KEY"],
    }],
    "temperature": 0.7,
    "timeout": 120,
}
```

### Azure OpenAI

```python
llm_config = {
    "config_list": [{
        "model": "gpt-4",
        "api_type": "azure",
        "api_base": "https://your-endpoint.openai.azure.com/",
        "api_key": os.environ["AZURE_OPENAI_KEY"],
        "api_version": "2023-05-15",
    }],
}
```

### 本地 LLM (LM Studio)

```python
llm_config = {
    "config_list": [{
        "model": "local-model",
        "base_url": "http://localhost:1234/v1",
        "api_key": "not-needed",
    }],
}
```

## 工作流示例

### 典型逆合成工作流

```
1. 用户提出目标化合物
   ↓
2. RetrosynthesisPlanner 分析目标
   - 调用 find_retrosynthesis_pathway()
   - 识别可能的前体
   ↓
3. EnzymeExpert 为每步选择酶
   - 调用 search_enzyme_by_ec()
   - 比较不同生物体的酶变体
   - 调用 compare_enzymes()
   ↓
4. KineticsAnalyst 评估效率
   - 调用 get_kinetic_parameters()
   - 计算反应速率
   - 识别瓶颈
   ↓
5. 团队综合结果
   - 推荐最优路径
   - 提供详细参数
   - 建议实验条件
```

## 性能优化

### 缓存知识库

```python
from agents import KnowledgeBaseTools

# 预加载知识库
tools = KnowledgeBaseTools("knowledge_base.jsonl")
tools._ensure_loaded()  # 一次性加载

# 多个智能体共享工具实例
agent1 = ChemoenzymaticAgent(kb_path="...", ...)
agent1.kb_tools = tools  # 共享

agent2 = ChemoenzymaticAgent(kb_path="...", ...)
agent2.kb_tools = tools  # 共享
```

### 限制上下文长度

```python
llm_config = {
    "config_list": [...],
    "max_tokens": 4000,  # 限制输出
}

# 在工具中限制返回数量
tools.search_enzyme_by_ec("1.1.1.1")  # 默认返回前 10 条
```

## 故障排除

### Autogen 未安装

```bash
pip install pyautogen
```

### API 密钥错误

```python
import os
print(os.environ.get("OPENAI_API_KEY"))  # 检查是否设置
```

### 知识库未找到

```bash
# 构建知识库
python -m knowledge_base.cli build --config kb_config.yaml
```

### Token 超限

- 使用 `gpt-3.5-turbo` 代替 `gpt-4`
- 减少工具返回的记录数量
- 限制对话轮数 (`max_round`)

## 进阶主题

### 集成 READRetro/RetroBioCat

```python
# 在智能体中添加专门工具调用 READRetro
def call_readretro(target_smiles: str) -> str:
    """Call READRetro for biosynthesis planning."""
    # 集成 READRetro 代码
    ...
    return results

agent.register_function(
    function_map={"call_readretro": call_readretro}
)
```

### 人类反馈循环

```python
from autogen import UserProxyAgent

user_proxy = UserProxyAgent(
    name="User",
    human_input_mode="ALWAYS",  # 每次都请求人类输入
    max_consecutive_auto_reply=0,
)

user_proxy.initiate_chat(agent, message="...")
# 智能体会在每步请求确认
```

### 成本控制

```python
# 使用更便宜的模型
llm_config = {
    "config_list": [
        {"model": "gpt-3.5-turbo", ...},  # 便宜
    ],
}

# 限制对话轮数
groupchat = GroupChat(
    agents=[...],
    max_round=10,  # 最多 10 轮
)

# 监控 token 使用
import autogen
autogen.logger.setLevel("INFO")  # 查看 token 统计
```

## 文档参考

- **Autogen 官方文档**: https://microsoft.github.io/autogen/
- **Knowledge Base API**: [../knowledge_base/README.md](../knowledge_base/README.md)
- **完整示例**: [../examples/autogen_retrosynthesis.py](../examples/autogen_retrosynthesis.py)

## License

遵循项目主 License
