#!/usr/bin/env python
"""
Production Retrosynthesis Agent - Main Entry Point

运行化酶逆合成规划Agent
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# 加载环境变量
load_dotenv()

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from agents.production_agent import RetrosynthesisAgent
from agents.logging_utils import setup_file_logging

logger = logging.getLogger(__name__)


def get_llm_config() -> dict:
    """
    从环境变量获取LLM配置
    支持: OpenAI, Azure OpenAI, OpenRouter
    """

    # 检查OpenRouter (推荐)
    if os.getenv("OPENROUTER_API_KEY"):
        model = os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5")
        logger.info(f"使用 OpenRouter API - 模型: {model}")
        return {
            "config_list": [
                {
                    "model": model,
                    "api_key": os.getenv("OPENROUTER_API_KEY"),
                    "base_url": "https://openrouter.ai/api/v1",
                }
            ],
            "temperature": float(os.getenv("LLM_TEMPERATURE", "0.7")),
            "max_tokens": int(os.getenv("LLM_MAX_TOKENS", "4000")),
        }

    # 检查OpenAI
    elif os.getenv("OPENAI_API_KEY"):
        logger.info("使用 OpenAI API")
        return {
            "config_list": [
                {
                    "model": os.getenv("OPENAI_MODEL", "gpt-4o-mini"),
                    "api_key": os.getenv("OPENAI_API_KEY"),
                }
            ],
            "temperature": float(os.getenv("LLM_TEMPERATURE", "0.7")),
            "max_tokens": int(os.getenv("LLM_MAX_TOKENS", "4000")),
        }

    # 检查Azure OpenAI
    elif os.getenv("AZURE_OPENAI_API_KEY"):
        logger.info("使用 Azure OpenAI API")
        return {
            "config_list": [
                {
                    "model": os.getenv("AZURE_OPENAI_DEPLOYMENT", "gpt-4"),
                    "api_type": "azure",
                    "api_key": os.getenv("AZURE_OPENAI_API_KEY"),
                    "base_url": os.getenv("AZURE_OPENAI_ENDPOINT"),
                    "api_version": os.getenv("AZURE_OPENAI_API_VERSION", "2024-02-01"),
                }
            ],
            "temperature": float(os.getenv("LLM_TEMPERATURE", "0.7")),
            "max_tokens": int(os.getenv("LLM_MAX_TOKENS", "4000")),
        }

    else:
        logger.error("未找到API密钥！请设置环境变量:")
        logger.error("  - OPENROUTER_API_KEY (推荐)")
        logger.error("  - OPENAI_API_KEY")
        logger.error("  - AZURE_OPENAI_API_KEY")
        sys.exit(1)


def interactive_mode(agent: RetrosynthesisAgent):
    """
    交互式模式
    """
    print("\n" + "="*60)
    print("🧬 化酶逆合成规划Agent - 交互模式")
    print("="*60)
    print("\n使用说明:")
    print("  - 输入目标分子名称或SMILES")
    print("  - 输入 'quick <分子>' 进行快速可行性评估")
    print("  - 输入 'quit' 退出")
    print("\n" + "="*60 + "\n")

    while True:
        try:
            user_input = input("\n请输入目标分子 (或命令): ").strip()

            if not user_input:
                continue

            if user_input.lower() in ['quit', 'exit', 'q']:
                print("\n再见! 👋")
                break

            # 快速检查模式
            if user_input.lower().startswith('quick '):
                target = user_input[6:].strip()
                print(f"\n正在快速评估 '{target}' ...\n")
                result = agent.quick_check(target)
                print("\n" + "="*60)
                print(result)
                print("="*60)

            # 完整规划模式
            else:
                target = user_input
                print(f"\n正在规划 '{target}' 的合成路线...\n")
                print("这可能需要1-2分钟，请稍候...\n")

                result = agent.plan(
                    target=target,
                    max_steps=6,
                    use_chemistry=True
                )

                print("\n" + "="*60)
                print("📋 规划结果")
                print("="*60 + "\n")
                print(result)
                print("\n" + "="*60 + "\n")

        except KeyboardInterrupt:
            print("\n\n收到中断信号，退出...")
            break
        except Exception as e:
            logger.error(f"错误: {e}")
            print(f"\n❌ 出错了: {e}")


def batch_mode(agent: RetrosynthesisAgent, targets: list, output_file: str):
    """
    批处理模式
    """
    print(f"\n批处理模式: 处理 {len(targets)} 个目标分子\n")

    results = []

    for i, target in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] 规划 '{target}' ...")

        try:
            result = agent.plan(target=target, max_steps=6, use_chemistry=True)
            results.append({
                "target": target,
                "status": "success",
                "result": result
            })
            print(f"  ✅ 完成\n")

        except Exception as e:
            logger.error(f"处理 '{target}' 时出错: {e}")
            results.append({
                "target": target,
                "status": "error",
                "error": str(e)
            })
            print(f"  ❌ 失败: {e}\n")

    # 保存结果
    import json
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n结果已保存到: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='化酶逆合成规划Agent',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 交互模式
  python run_agent.py

  # 规划单个分子
  python run_agent.py --target "布洛芬"

  # 快速检查
  python run_agent.py --target "L-DOPA" --quick

  # 批处理
  python run_agent.py --batch targets.txt --output results.json

  # 指定知识库
  python run_agent.py --kb my_kb.jsonl --target "阿司匹林"

  # 启用详细日志
  python run_agent.py --target "布洛芬" --debug

  # 自定义日志目录
  python run_agent.py --target "布洛芬" --log-dir my_logs
        """
    )

    parser.add_argument(
        '--kb',
        type=str,
        default='knowledge_base_output/knowledge_base.jsonl',
        help='知识库文件路径 (默认: knowledge_base_output/knowledge_base.jsonl)'
    )

    parser.add_argument(
        '--target',
        type=str,
        help='目标分子 (名称或SMILES)'
    )

    parser.add_argument(
        '--quick',
        action='store_true',
        help='快速可行性检查 (不做完整规划)'
    )

    parser.add_argument(
        '--batch',
        type=str,
        help='批处理模式: 从文件读取目标列表 (每行一个)'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='batch_results.json',
        help='批处理结果输出文件 (默认: batch_results.json)'
    )

    parser.add_argument(
        '--max-steps',
        type=int,
        default=6,
        help='最大合成步数 (默认: 6)'
    )

    parser.add_argument(
        '--no-chemistry',
        action='store_true',
        help='禁用传统化学反应 (仅生物催化)'
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help='启用DEBUG级别日志'
    )

    parser.add_argument(
        '--log-dir',
        type=str,
        default='logs',
        help='日志目录 (默认: logs)'
    )

    parser.add_argument(
        '--no-session-log',
        action='store_true',
        help='禁用会话日志 (工具调用追踪)'
    )

    args = parser.parse_args()

    # 配置日志级别和文件输出
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 降低第三方库的日志级别，减少噪音
    logging.getLogger('httpx').setLevel(logging.WARNING)
    logging.getLogger('httpcore').setLevel(logging.WARNING)
    logging.getLogger('openai').setLevel(logging.WARNING)
    logging.getLogger('h5py').setLevel(logging.WARNING)
    logging.getLogger('tensorflow').setLevel(logging.WARNING)
    logging.getLogger('pubchempy').setLevel(logging.INFO)  # 避免DEBUG信息过多

    # 设置autogen日志级别 - 避免重复的成本警告
    logging.getLogger('autogen').setLevel(logging.WARNING)
    logging.getLogger('autogen.oai.client').setLevel(logging.ERROR)  # 屏蔽模型成本警告

    # 设置文件日志
    log_dir = Path(args.log_dir)
    log_file = setup_file_logging(log_dir, level=log_level)
    logger.info(f"日志文件: {log_file}")

    # 检查知识库
    kb_path = Path(args.kb)
    if not kb_path.exists():
        logger.error(f"知识库文件不存在: {kb_path}")
        logger.error("请先构建知识库: uv run python -m knowledge_base.cli build")
        sys.exit(1)

    # 获取LLM配置
    llm_config = get_llm_config()

    # 创建Agent (带日志配置)
    logger.info("正在初始化Agent...")
    enable_session_log = not args.no_session_log
    session_log_dir = log_dir / "sessions"

    agent = RetrosynthesisAgent(
        kb_path=kb_path,
        llm_config=llm_config,
        name="化酶逆合成规划师",
        enable_logging=enable_session_log,
        log_dir=session_log_dir
    )
    logger.info("Agent初始化完成")

    if enable_session_log:
        logger.info(f"会话日志目录: {session_log_dir}")
        logger.info(f"当前会话日志: {agent.session_logger.get_session_file()}")

    # 选择运行模式
    if args.batch:
        # 批处理模式
        batch_file = Path(args.batch)
        if not batch_file.exists():
            logger.error(f"批处理文件不存在: {batch_file}")
            sys.exit(1)

        with open(batch_file) as f:
            targets = [line.strip() for line in f if line.strip()]

        batch_mode(agent, targets, args.output)

    elif args.target:
        # 单分子模式
        if args.quick:
            print(f"\n正在快速评估 '{args.target}' ...\n")
            result = agent.quick_check(args.target)
        else:
            print(f"\n正在规划 '{args.target}' 的合成路线...\n")
            result = agent.plan(
                target=args.target,
                max_steps=args.max_steps,
                use_chemistry=not args.no_chemistry
            )

        print("\n" + "="*60)
        print(result)
        print("="*60 + "\n")

    else:
        # 交互模式
        interactive_mode(agent)


if __name__ == "__main__":
    main()
