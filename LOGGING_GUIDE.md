# Agent 日志系统使用指南

本项目现已集成全面的日志系统，方便调试和追踪 Agent 的工作流程。

## 日志类型

### 1. 标准日志 (logs/*.log)
记录程序运行的基本信息、错误等。

**位置**: `logs/agent_YYYYMMDD_HHMMSS.log`

**内容**:
- 程序启动/初始化信息
- 工具调用的简要提示
- 错误和警告信息
- 系统状态

### 2. 会话日志 (logs/sessions/*.jsonl)
详细记录每次运行的所有工具调用、LLM消息、统计信息。

**位置**: `logs/sessions/session_YYYYMMDD_HHMMSS.jsonl`

**内容** (JSONL格式，每行一个JSON对象):
- `tool_call`: 工具调用详情（输入、输出、耗时、错误）
- `llm_message`: LLM对话消息
- `session_summary`: 会话摘要和统计

## 使用方法

### 基本用法

```bash
# 默认启用日志
uv run python run_agent.py --target "布洛芬"

# 查看日志位置
# 标准日志: logs/agent_20251210_142030.log
# 会话日志: logs/sessions/session_20251210_142030.jsonl
```

### 调试模式

```bash
# 启用DEBUG级别日志 (更详细)
uv run python run_agent.py --target "布洛芬" --debug
```

### 自定义日志目录

```bash
# 指定日志目录
uv run python run_agent.py --target "布洛芬" --log-dir my_logs
```

### 禁用会话日志

```bash
# 只保留标准日志，不记录工具调用详情
uv run python run_agent.py --target "布洛芬" --no-session-log
```

## 查看日志

### 方法1: 直接查看文件

```bash
# 查看标准日志
tail -f logs/agent_*.log

# 查看会话日志 (JSONL格式)
cat logs/sessions/session_*.jsonl | jq .
```

### 方法2: 使用分析工具

```python
from agents.logging_utils import parse_session_log, print_session_analysis

# 解析会话日志
session_file = "logs/sessions/session_20251210_142030.jsonl"
data = parse_session_log(session_file)

# 打印分析报告
print_session_analysis(session_file)
```

输出示例:
```
======================================================================
📊 Agent会话分析报告
======================================================================

目标分子: ibuprofen
状态: ✅ 成功

总工具调用: 4
总LLM调用: 2
总耗时: 13.120s

工具使用统计:
  - name_to_smiles: 1次
  - plan_biocatalytic_route: 1次
  - search_enzyme_by_ec: 1次
  - get_kinetic_parameters: 1次

工具调用详情 (共4次):
  [1] ✅ name_to_smiles - 0.15s
  [2] ✅ plan_biocatalytic_route - 12.3s
  [3] ✅ search_enzyme_by_ec - 0.42s
  [4] ❌ get_kinetic_parameters - 0.25s
      错误: EC number not found in database
======================================================================
```

## 会话日志格式

### 工具调用记录
```json
{
  "timestamp": "2025-12-10T14:28:14.952840",
  "type": "tool_call",
  "tool_name": "name_to_smiles",
  "inputs": {
    "args": ["ibuprofen"],
    "kwargs": {}
  },
  "outputs": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "duration": 0.15,
  "error": null,
  "success": true
}
```

### LLM消息记录
```json
{
  "timestamp": "2025-12-10T14:28:14.959017",
  "type": "llm_message",
  "role": "user",
  "content": "请帮我设计布洛芬的化酶合成路线。",
  "name": null,
  "function_call": null
}
```

### 会话摘要
```json
{
  "timestamp": "2025-12-10T14:28:14.962095",
  "type": "session_summary",
  "target": "ibuprofen",
  "result": "路线规划完成，共找到2条可行路线",
  "success": true,
  "statistics": {
    "total_tool_calls": 4,
    "total_llm_calls": 2,
    "total_time": 13.12,
    "tool_breakdown": {
      "name_to_smiles": 1,
      "plan_biocatalytic_route": 1,
      "search_enzyme_by_ec": 1,
      "get_kinetic_parameters": 1
    }
  }
}
```

## 调试技巧

### 1. 追踪工具调用流程

查看 Agent 调用了哪些工具、顺序如何：

```bash
cat logs/sessions/session_*.jsonl | jq -r 'select(.type=="tool_call") | "\(.tool_name) - \(.duration)s"'
```

输出:
```
name_to_smiles - 0.15s
plan_biocatalytic_route - 12.3s
search_enzyme_by_ec - 0.42s
get_kinetic_parameters - 0.25s
```

### 2. 查找错误

找出失败的工具调用：

```bash
cat logs/sessions/session_*.jsonl | jq 'select(.type=="tool_call" and .success==false)'
```

### 3. 统计工具使用

```bash
cat logs/sessions/session_*.jsonl | \
  jq -r 'select(.type=="tool_call") | .tool_name' | \
  sort | uniq -c | sort -rn
```

### 4. 查看具体工具的输入输出

```bash
# 查看 plan_biocatalytic_route 的输入和输出
cat logs/sessions/session_*.jsonl | \
  jq 'select(.tool_name=="plan_biocatalytic_route") | {inputs, outputs}'
```

## 编程接口

### 在代码中使用会话日志

```python
from agents.production_agent import RetrosynthesisAgent

# 创建 Agent (默认启用日志)
agent = RetrosynthesisAgent(
    kb_path="knowledge_base_output/knowledge_base.jsonl",
    llm_config=llm_config,
    enable_logging=True,  # 启用日志
    log_dir="logs/sessions"  # 日志目录
)

# 运行规划
result = agent.plan(target="ibuprofen")

# 获取当前会话日志文件
if agent.session_logger:
    log_file = agent.session_logger.get_session_file()
    print(f"会话日志: {log_file}")
```

### 自定义工具日志

```python
from agents.logging_utils import AgentSessionLogger

# 创建日志记录器
logger = AgentSessionLogger(session_dir="my_logs")

# 手动记录工具调用
logger.log_tool_call(
    tool_name="my_custom_tool",
    inputs={"param1": "value1"},
    outputs={"result": "data"},
    duration=1.23,
    error=None
)

# 记录LLM消息
logger.log_llm_message(
    role="assistant",
    content="Response from LLM"
)

# 记录会话摘要
logger.log_session_summary(
    target="molecule_name",
    result="Success",
    success=True
)
```

## 测试日志系统

运行测试脚本验证日志功能：

```bash
uv run python test_logging.py
```

这会创建一个测试会话，记录模拟的工具调用和LLM消息，并展示分析报告。

## 常见问题

### Q: 日志文件太多怎么办？

可以定期清理旧日志：

```bash
# 删除7天前的日志
find logs -name "*.log" -mtime +7 -delete
find logs/sessions -name "*.jsonl" -mtime +7 -delete
```

### Q: 如何只看特定工具的调用？

```bash
cat logs/sessions/session_*.jsonl | \
  jq 'select(.tool_name=="search_enzyme_by_ec")'
```

### Q: 日志文件在哪里？

- 标准日志: `logs/agent_YYYYMMDD_HHMMSS.log`
- 会话日志: `logs/sessions/session_YYYYMMDD_HHMMSS.jsonl`

运行时会在控制台打印日志文件位置。

## 最佳实践

1. **调试时开启DEBUG模式**: `--debug` 获取更详细的日志
2. **性能分析**: 查看会话日志中各工具的 `duration` 找出瓶颈
3. **错误追踪**: 搜索 `"success": false` 快速定位问题
4. **定期备份**: 重要的会话日志应该备份保存
5. **使用 jq 处理**: `jq` 是处理 JSONL 格式的利器

## 工具列表

所有工具调用都会被记录：

### 分子处理
- `name_to_smiles`
- `analyze_molecule`

### 路线规划 (RetroBioCat)
- `plan_biocatalytic_route`
- `find_enzymatic_reactions`
- `check_commercial_availability`
- `compare_retrosynthesis_approaches`
- `analyze_pathway_feasibility`

### 知识库查询
- `search_enzyme_by_ec`
- `search_reactions_by_compound`
- `find_retrosynthesis_pathway`
- `get_kinetic_parameters`
- `get_enzyme_sequence`
- `compare_enzymes`
- `get_statistics`
