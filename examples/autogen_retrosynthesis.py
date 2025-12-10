#!/usr/bin/env python3
"""
Example: Using Autogen AI Agents for Chemoenzymatic Retrosynthesis

This script demonstrates:
1. Single-agent retrosynthesis planning with RAG
2. Multi-agent team collaboration
3. Enzyme selection and optimization
4. Custom retrosynthesis queries

Setup:
1. Copy .env.example to .env
2. Fill in your API keys in .env:
   - For OpenRouter: OPENROUTER_API_KEY, OPENROUTER_BASE_URL, OPENROUTER_BASE_MODEL
   - Or for OpenAI: OPENAI_API_KEY
3. Run: python examples/autogen_retrosynthesis.py
"""

import os
from pathlib import Path

from agents import (
    ChemoenzymaticAgent,
    create_retrosynthesis_team,
    get_default_config,
    get_openrouter_config,
    print_config_status,
)


def example_1_single_agent_planning():
    """Example 1: Single agent retrosynthesis planning."""
    print("=" * 70)
    print("Example 1: Single Agent Retrosynthesis Planning")
    print("=" * 70)

    # Configuration - automatically loads from .env
    kb_path = "knowledge_base_output/knowledge_base.jsonl"
    llm_config = get_default_config(temperature=0.7)

    print(f"\nUsing model: {llm_config['config_list'][0]['model']}")

    # Create agent
    agent = ChemoenzymaticAgent(
        kb_path=kb_path,
        llm_config=llm_config,
        name="RetrosynthesisExpert",
    )

    # Query: Design synthesis route
    query = """Design an enzymatic synthesis route for glucose-6-phosphate.

Available starting materials: glucose, ATP

Requirements:
- Prefer E. coli enzymes
- Provide kinetic parameters
- Suggest alternative routes if available
"""

    print(f"\nQuery: {query}")
    print("\nAgent working...\n")

    # Chat with agent
    response = agent.chat(query)
    print(f"\nResponse: {response}\n")


def example_2_multi_agent_team():
    """Example 2: Multi-agent team for complex retrosynthesis."""
    print("\n" + "=" * 70)
    print("Example 2: Multi-Agent Team Collaboration")
    print("=" * 70)

    kb_path = "knowledge_base_output/knowledge_base.jsonl"
    llm_config = get_default_config(temperature=0.7)

    print(f"\nUsing model: {llm_config['config_list'][0]['model']}")

    # Create multi-agent team
    team = create_retrosynthesis_team(kb_path, llm_config)

    # Complex query requiring collaboration
    query = """I need to synthesize L-DOPA (L-3,4-dihydroxyphenylalanine) using enzymatic methods.

Target: L-DOPA
Available substrates: tyrosine, phenylalanine, glucose

Requirements:
1. Design a multi-step enzymatic pathway
2. Compare enzyme options for each step
3. Provide kinetic analysis
4. Evaluate overall pathway efficiency
5. Suggest organism engineering strategy

The team should collaborate:
- RetrosynthesisPlanner: Design the pathway
- EnzymeExpert: Select optimal enzymes
- KineticsAnalyst: Analyze reaction rates
"""

    print(f"\nComplex Query: {query}")
    print("\nTeam collaborating...\n")

    # Initiate group chat
    user_proxy = team["user_proxy"]
    manager = team["manager"]

    user_proxy.initiate_chat(manager, message=query)

    print("\nTeam discussion completed!\n")


def example_3_high_level_session():
    """Example 3: Using high-level RetrosynthesisSession."""
    print("\n" + "=" * 70)
    print("Example 3: High-Level Retrosynthesis Session")
    print("=" * 70)

    kb_path = "knowledge_base_output/knowledge_base.jsonl"
    llm_config = get_default_config(temperature=0.7)

    print(f"\nUsing model: {llm_config['config_list'][0]['model']}")

    # Create session with multi-agent team
    session = RetrosynthesisSession(
        kb_path=kb_path,
        llm_config=llm_config,
        use_multi_agent=True,
    )

    # Example 3a: Plan synthesis
    print("\n--- Task 1: Plan Synthesis ---")
    session.plan_synthesis(
        target_compound="acetaldehyde",
        available_substrates=["ethanol", "NAD+"],
        constraints={
            "organism": "Saccharomyces cerevisiae",
            "temperature_max": 37,
        },
    )

    # Example 3b: Select enzyme
    print("\n--- Task 2: Select Optimal Enzyme ---")
    session.select_enzyme(
        ec_number="1.1.1.1",  # Alcohol dehydrogenase
        criteria={
            "kcat_min": 100,  # Minimum turnover rate
            "organism_preference": "E. coli",
            "has_sequence": True,
        },
    )

    # Example 3c: Optimize conditions
    print("\n--- Task 3: Optimize Reaction Conditions ---")
    session.optimize_conditions(
        ec_number="1.1.1.1",
        substrate="ethanol",
    )


def example_4_direct_tool_usage():
    """Example 4: Direct use of KB tools without agents."""
    print("\n" + "=" * 70)
    print("Example 4: Direct Knowledge Base Tool Usage")
    print("=" * 70)

    from agents import KnowledgeBaseTools

    kb_path = "knowledge_base_output/knowledge_base.jsonl"
    tools = KnowledgeBaseTools(kb_path)

    # Search enzyme by EC
    print("\n1. Search enzyme EC 1.1.1.1:")
    result = tools.search_enzyme_by_ec("1.1.1.1")
    print(result[:500] + "...")

    # Find retrosynthesis pathway
    print("\n2. Find pathway for glucose-6-phosphate:")
    result = tools.find_retrosynthesis_pathway(
        target_compound="glucose-6-phosphate",
        available_substrates=["glucose", "ATP"],
    )
    print(result[:500] + "...")

    # Get kinetic parameters
    print("\n3. Get kinetics for alcohol dehydrogenase:")
    result = tools.get_kinetic_parameters(
        ec_number="1.1.1.1",
        organism="Escherichia coli",
    )
    print(result[:500] + "...")

    # Compare enzymes
    print("\n4. Compare multiple alcohol dehydrogenases:")
    result = tools.compare_enzymes(["1.1.1.1", "1.1.1.2"])
    print(result[:500] + "...")


def example_5_custom_agent():
    """Example 5: Create custom agent with specific expertise."""
    print("\n" + "=" * 70)
    print("Example 5: Custom Agent - Protein Engineering Expert")
    print("=" * 70)

    from autogen import AssistantAgent, UserProxyAgent
    from agents import KnowledgeBaseTools

    kb_path = "knowledge_base_output/knowledge_base.jsonl"
    kb_tools = KnowledgeBaseTools(kb_path)

    llm_config = get_default_config(temperature=0.7)

    print(f"\nUsing model: {llm_config['config_list'][0]['model']}")

    # Create specialized agent
    protein_engineer = AssistantAgent(
        name="ProteinEngineeringExpert",
        system_message="""You are a protein engineering expert specializing in enzyme optimization.

Your expertise:
- Analyzing enzyme sequences for mutation sites
- Predicting structure-function relationships
- Designing directed evolution strategies
- Evaluating enzyme promiscuity

Use the knowledge base to retrieve enzyme sequences and compare variants.""",
        llm_config=llm_config,
    )

    # Register tools
    protein_engineer.register_function(
        function_map={
            "get_enzyme_sequence": kb_tools.get_enzyme_sequence,
            "search_enzyme_by_ec": kb_tools.search_enzyme_by_ec,
            "compare_enzymes": kb_tools.compare_enzymes,
        }
    )

    # Create user proxy
    user = UserProxyAgent(
        name="User",
        human_input_mode="NEVER",
        max_consecutive_auto_reply=5,
        code_execution_config=False,
    )

    # Query
    query = """I want to engineer an alcohol dehydrogenase (EC 1.1.1.1) for improved activity on bulky substrates.

Tasks:
1. Retrieve sequences from different organisms
2. Compare sequence conservation
3. Identify potential mutation sites
4. Suggest a directed evolution strategy
"""

    print(f"\nQuery: {query}")
    print("\nProtein engineer working...\n")

    user.initiate_chat(protein_engineer, message=query)


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("Autogen Chemoenzymatic Retrosynthesis Examples")
    print("=" * 70)

    # Print configuration status
    print_config_status()

    # Check if any LLM API is configured
    try:
        config = get_default_config()
        print(f"\n✓ Using LLM: {config['config_list'][0]['model']}\n")
    except ValueError as e:
        print(f"\n⚠️  Warning: {e}")
        print("\nRunning example 4 (no API key required)...\n")
        example_4_direct_tool_usage()
        return

    # Check KB exists
    kb_path = Path("knowledge_base_output/knowledge_base.jsonl")
    if not kb_path.exists():
        print(f"\n⚠️  Knowledge base not found: {kb_path}")
        print("Please build it first:")
        print("  python -m knowledge_base.cli build --config kb_config.yaml")
        return

    # Run examples
    try:
        # Example 4 doesn't need LLM
        example_4_direct_tool_usage()

        # Examples with LLM
        example_1_single_agent_planning()

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
