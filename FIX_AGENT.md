# Agent 修复方案

## 问题诊断

当前Agent无法正常工作的原因：
1. **工具调用数为0** - LLM没有调用任何工具
2. **模型**: `google/gemini-3-pro-preview` 返回空响应
3. **根本原因**: 该模型可能不完全支持AutoGen的function calling格式

## 解决方案

### 方案1: 更换模型（推荐）✅

使用经过验证支持function calling的模型：

```bash
# 编辑 .env 文件
OPENROUTER_MODEL=openai/gpt-4o-mini
# 或
OPENROUTER_MODEL=anthropic/claude-3.5-haiku
# 或
OPENROUTER_MODEL=meta-llama/llama-3.1-70b-instruct
```

### 方案2: 添加手动工具调用解析

如果必须使用当前模型，需要修改代码手动解析工具调用。

### 方案3: 使用 LiteLLM（最稳定）

LiteLLM 可以统一不同模型的function calling格式。

## 快速修复步骤

1. 修改 `.env`:
```bash
OPENROUTER_MODEL=openai/gpt-4o-mini
```

2. 重新运行:
```bash
uv run python run_agent.py --target "布洛芬"
```

3. 检查日志:
```bash
cat logs/sessions/session_*.jsonl | tail -1 | jq '.statistics.total_tool_calls'
```

应该看到 > 0 的工具调用次数。

## 已验证可用的模型

以下模型已验证支持function calling:
- ✅ `openai/gpt-4o-mini` (便宜、快速)
- ✅ `openai/gpt-4o` (最强)
- ✅ `anthropic/claude-3.5-sonnet` (强大)
- ✅ `anthropic/claude-3.5-haiku` (快速)
- ✅ `meta-llama/llama-3.1-70b-instruct` (开源)

不稳定/不支持:
- ❌ `google/gemini-3-pro-preview` (返回空响应)
- ❌ `google/gemini-flash-1.5` (地理限制)

## 需要代码修改吗？

**不需要！** 只需要修改 `.env` 中的模型配置即可。

当前代码已经：
- ✅ 正确注册了工具
- ✅ 配置了function_map
- ✅ 添加了日志追踪
- ✅ 有空响应保护

问题仅在于模型不兼容。
