#!/usr/bin/env python
"""
RetroBioCat 工具基本使用示例

展示如何使用配置好的 RetroBioCat 工具进行逆合成规划

运行方式：
    从项目根目录运行: uv run python examples/basic_usage.py
"""

import json
import os
import sys

# 添加项目根目录到路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

os.environ['TOKENIZERS_PARALLELISM'] = 'false'

from agents.retrobiocat_tools import RetroBioCatTools
from agents.utils import name_to_smiles


def example_1_name_to_smiles():
    """示例1: 分子名称转 SMILES"""
    print("=" * 60)
    print("示例 1: 分子名称转 SMILES")
    print("=" * 60)

    molecules = ["aspirin", "ibuprofen", "caffeine", "glucose"]

    for name in molecules:
        result = json.loads(name_to_smiles(name))
        if result.get("success"):
            print(f"  {name}: {result['smiles']}")
        else:
            print(f"  {name}: 转换失败 - {result.get('error')}")
    print()


def example_2_find_reactions():
    """示例2: 查找单步酶促反应"""
    print("=" * 60)
    print("示例 2: 查找单步酶促反应")
    print("=" * 60)

    rbc = RetroBioCatTools()

    # 测试乙醇
    target = "CCO"  # 乙醇
    print(f"目标分子: 乙醇 (SMILES: {target})")

    result = json.loads(rbc.find_enzymatic_reactions(target))

    print(f"找到 {result['total_reactions']} 个反应")
    print(f"各数据库结果: {result['by_expander']}")
    print()

    # 显示前 3 个反应
    print("Top 3 反应:")
    for rxn in result["reactions"][:3]:
        print(f"  #{rxn['rank']}: {rxn['name']}")
        print(f"      类型: {rxn['domain']} - {rxn['type']}")
        print(f"      评分: {rxn['score']}")
        if rxn.get("precedents"):
            print(f"      有 {len(rxn['precedents'])} 个文献先例")
        print()


def example_3_check_availability():
    """示例3: 检查商业可获得性"""
    print("=" * 60)
    print("示例 3: 检查商业可获得性")
    print("=" * 60)

    rbc = RetroBioCatTools()

    molecules = [
        ("CCO", "乙醇"),
        ("CC(=O)O", "乙酸"),
        ("c1ccccc1", "苯"),
        ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "布洛芬"),
    ]

    smiles_list = [m[0] for m in molecules]
    result = json.loads(rbc.check_commercial_availability(smiles_list))

    print(f"检查 {result['total_molecules']} 个分子")
    print(f"可购买: {result['available_count']} 个")
    print()

    for i, mol_data in enumerate(result["molecules"]):
        name = molecules[i][1]
        status = "✅ 可购买" if mol_data["available"] else "❌ 不可购买"
        print(f"  {name}: {status}")
    print()


def example_4_plan_route():
    """示例4: 规划混合逆合成路线"""
    print("=" * 60)
    print("示例 4: 规划混合逆合成路线 (MCTS)")
    print("=" * 60)

    rbc = RetroBioCatTools()

    # 简单目标分子
    target = "CCO"  # 乙醇
    print(f"目标分子: 乙醇 (SMILES: {target})")
    print("配置: max_steps=3, max_search_time=20s, 混合模式")
    print("搜索中...")

    result = json.loads(rbc.plan_biocatalytic_route(
        target_smiles=target,
        max_steps=3,
        use_chemistry=True,
        max_search_time=20
    ))

    print(f"找到 {result['pathways_found']} 条路径")
    print()

    if result["pathways_found"] > 0:
        pathway = result["pathways"][0]
        print("最佳路径:")
        print(f"  总步数: {pathway['total_steps']}")
        print(f"  生物催化步骤: {pathway['biocatalysis_steps']}")
        print(f"  化学步骤: {pathway['chemistry_steps']}")
        print()

        print("  反应序列:")
        for rxn in pathway["reactions"]:
            print(f"    Step {rxn['step']}: {rxn['name']}")
            print(f"      类型: {rxn['domain']}")
            print(f"      评分: {rxn['score']}")
        print()

        print("  起始原料:")
        for sm in pathway["starting_materials"]:
            status = "✅" if sm["commercially_available"] else "❌"
            print(f"    {status} {sm['smiles']}")
    print()


def example_5_bio_only():
    """示例5: 仅生物催化规划"""
    print("=" * 60)
    print("示例 5: 仅生物催化规划")
    print("=" * 60)

    rbc = RetroBioCatTools()

    target = "CC(C)CC(C(=O)O)N"  # 亮氨酸
    print(f"目标分子: 亮氨酸")
    print(f"SMILES: {target}")
    print("配置: max_steps=4, 仅生物催化")
    print("搜索中...")

    result = json.loads(rbc.plan_biocatalytic_route(
        target_smiles=target,
        max_steps=4,
        use_chemistry=False,  # 仅生物催化
        max_search_time=30
    ))

    print(f"找到 {result['pathways_found']} 条路径")

    if result["pathways_found"] > 0:
        pathway = result["pathways"][0]
        print(f"最佳路径: {pathway['total_steps']} 步")
        print(f"全部使用生物催化: {pathway['biocatalysis_steps']} 步")
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("RetroBioCat 工具使用示例")
    print("=" * 60 + "\n")

    example_1_name_to_smiles()
    example_2_find_reactions()
    example_3_check_availability()
    example_4_plan_route()
    example_5_bio_only()

    print("=" * 60)
    print("✅ 所有示例运行完成!")
    print("=" * 60)
