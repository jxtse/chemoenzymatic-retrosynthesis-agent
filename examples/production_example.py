#!/usr/bin/env python
"""
Production Agent Example - 生产环境Agent使用示例

展示如何使用RetrosynthesisAgent进行逆合成规划
"""

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# 加载环境变量
load_dotenv()

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from agents.production_agent import RetrosynthesisAgent


def example_1_basic_planning():
    """
    示例1: 基本的逆合成规划
    """
    print("\n" + "="*60)
    print("示例1: 基本逆合成规划 - 布洛芬")
    print("="*60 + "\n")

    # LLM配置
    llm_config = {
        "config_list": [
            {
                "model": os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5"),
                "api_key": os.getenv("OPENROUTER_API_KEY"),
                "base_url": "https://openrouter.ai/api/v1",
            }
        ],
        "temperature": 0.7,
        "max_tokens": 4000,
    }

    # 创建Agent
    agent = RetrosynthesisAgent(
        kb_path="knowledge_base_output/knowledge_base.jsonl",
        llm_config=llm_config,
    )

    # 规划路线
    result = agent.plan(
        target="ibuprofen",  # 布洛芬
        max_steps=6,
        use_chemistry=True
    )

    print(result)


def example_2_quick_check():
    """
    示例2: 快速可行性检查
    """
    print("\n" + "="*60)
    print("示例2: 快速可行性检查 - L-DOPA")
    print("="*60 + "\n")

    llm_config = {
        "config_list": [
            {
                "model": os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5"),
                "api_key": os.getenv("OPENROUTER_API_KEY"),
                "base_url": "https://openrouter.ai/api/v1",
            }
        ],
        "temperature": 0.7,
    }

    agent = RetrosynthesisAgent(
        kb_path="knowledge_base_output/knowledge_base.jsonl",
        llm_config=llm_config,
    )

    # 快速检查
    result = agent.quick_check("L-DOPA")
    print(result)


def example_3_bio_only():
    """
    示例3: 纯生物催化路线
    """
    print("\n" + "="*60)
    print("示例3: 纯生物催化路线 - 苯乙胺")
    print("="*60 + "\n")

    llm_config = {
        "config_list": [
            {
                "model": os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5"),
                "api_key": os.getenv("OPENROUTER_API_KEY"),
                "base_url": "https://openrouter.ai/api/v1",
            }
        ],
        "temperature": 0.7,
    }

    agent = RetrosynthesisAgent(
        kb_path="knowledge_base_output/knowledge_base.jsonl",
        llm_config=llm_config,
    )

    # 禁用化学反应，仅使用生物催化
    result = agent.plan(
        target="phenethylamine",  # 苯乙胺
        max_steps=5,
        use_chemistry=False  # 纯生物催化
    )

    print(result)


def example_4_batch_processing():
    """
    示例4: 批量处理多个目标
    """
    print("\n" + "="*60)
    print("示例4: 批量处理")
    print("="*60 + "\n")

    llm_config = {
        "config_list": [
            {
                "model": os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5"),
                "api_key": os.getenv("OPENROUTER_API_KEY"),
                "base_url": "https://openrouter.ai/api/v1",
            }
        ],
        "temperature": 0.7,
    }

    agent = RetrosynthesisAgent(
        kb_path="knowledge_base_output/knowledge_base.jsonl",
        llm_config=llm_config,
    )

    # 批量目标
    targets = [
        "aspirin",      # 阿司匹林
        "paracetamol",  # 扑热息痛
        "caffeine",     # 咖啡因
    ]

    results = []

    for target in targets:
        print(f"\n处理: {target}")
        try:
            result = agent.quick_check(target)
            results.append({
                "target": target,
                "status": "success",
                "result": result
            })
            print(f"✅ 完成")
        except Exception as e:
            print(f"❌ 失败: {e}")
            results.append({
                "target": target,
                "status": "error",
                "error": str(e)
            })

    # 保存结果
    import json
    with open("batch_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n结果已保存到 batch_results.json")


def example_5_from_config():
    """
    示例5: 从配置文件创建Agent
    """
    print("\n" + "="*60)
    print("示例5: 从配置文件创建")
    print("="*60 + "\n")

    from agents.production_agent import create_agent_from_config

    # 创建临时配置
    import json
    config = {
        "name": "我的逆合成Agent",
        "model": os.getenv("OPENROUTER_MODEL", "google/gemini-flash-1.5"),
        "api_key": os.getenv("OPENROUTER_API_KEY"),
        "temperature": 0.7,
        "max_tokens": 4000,
    }

    with open("temp_config.json", "w") as f:
        json.dump(config, f, indent=2)

    # 从配置创建
    agent = create_agent_from_config(
        config_path="temp_config.json",
        kb_path="knowledge_base_output/knowledge_base.jsonl"
    )

    result = agent.quick_check("vanillin")  # 香草醛
    print(result)

    # 清理
    os.remove("temp_config.json")


def main():
    """
    主函数: 运行所有示例
    """
    import argparse

    parser = argparse.ArgumentParser(description="Production Agent示例")
    parser.add_argument(
        "--example",
        type=int,
        choices=[1, 2, 3, 4, 5],
        help="运行指定示例 (1-5)"
    )

    args = parser.parse_args()

    examples = {
        1: ("基本逆合成规划", example_1_basic_planning),
        2: ("快速可行性检查", example_2_quick_check),
        3: ("纯生物催化路线", example_3_bio_only),
        4: ("批量处理", example_4_batch_processing),
        5: ("从配置文件创建", example_5_from_config),
    }

    if args.example:
        # 运行指定示例
        name, func = examples[args.example]
        print(f"\n运行示例 {args.example}: {name}")
        func()
    else:
        # 运行所有示例
        print("\n" + "="*60)
        print("Production Agent 示例集")
        print("="*60)
        print("\n可用示例:")
        for num, (name, _) in examples.items():
            print(f"  {num}. {name}")
        print("\n使用 --example N 运行指定示例")
        print("例如: python production_example.py --example 1")


if __name__ == "__main__":
    main()
