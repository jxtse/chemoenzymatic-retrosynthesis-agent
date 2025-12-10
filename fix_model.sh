#!/bin/bash
# 自动修复Agent模型配置

echo "🔧 Agent 模型修复脚本"
echo "================================"
echo ""

# 检查当前模型
CURRENT_MODEL=$(grep "^OPENROUTER_MODEL=" .env 2>/dev/null | cut -d'=' -f2)
echo "当前模型: $CURRENT_MODEL"
echo ""

# 推荐模型列表
echo "推荐的可用模型："
echo "  1) openai/gpt-4o-mini (推荐: 便宜、快速、可靠)"
echo "  2) anthropic/claude-3.5-haiku (快速)"
echo "  3) anthropic/claude-3.5-sonnet (强大)"
echo "  4) meta-llama/llama-3.1-70b-instruct (开源)"
echo "  5) 保持当前模型"
echo ""

read -p "选择模型 [1-5]: " choice

case $choice in
  1)
    NEW_MODEL="openai/gpt-4o-mini"
    ;;
  2)
    NEW_MODEL="anthropic/claude-3.5-haiku"
    ;;
  3)
    NEW_MODEL="anthropic/claude-3.5-sonnet"
    ;;
  4)
    NEW_MODEL="meta-llama/llama-3.1-70b-instruct"
    ;;
  5)
    echo "保持当前配置"
    exit 0
    ;;
  *)
    echo "无效选择"
    exit 1
    ;;
esac

echo ""
echo "正在更新模型为: $NEW_MODEL"

# 备份.env
cp .env .env.backup.$(date +%Y%m%d_%H%M%S)

# 更新模型
if grep -q "^OPENROUTER_MODEL=" .env; then
  sed -i "s|^OPENROUTER_MODEL=.*|OPENROUTER_MODEL=$NEW_MODEL|" .env
else
  echo "OPENROUTER_MODEL=$NEW_MODEL" >> .env
fi

echo "✅ 模型已更新!"
echo ""
echo "测试新配置:"
echo "  uv run python run_agent.py --target '布洛芬'"
echo ""
