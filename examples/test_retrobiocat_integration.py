#!/usr/bin/env python3
"""
Test RetroBioCat Integration with Chemoenzymatic Agent

This script tests the integration of RetroBioCat 2.0 tools into the agent.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_retrobiocat_tools_standalone():
    """Test RetroBioCat tools without agent."""
    print("=" * 70)
    print("Test 1: RetroBioCat Tools (Standalone)")
    print("=" * 70)

    from agents import RetroBioCatTools

    tools = RetroBioCatTools()

    # Test 1: Find enzymatic reactions
    print("\n--- Finding enzymatic reactions for phenylethanol ---")
    target = "c1ccc(cc1)CCO"  # Phenylethanol

    try:
        result = tools.find_enzymatic_reactions(
            target,
            expander_types=['retrobiocat']
        )
        print(result[:500] + "..." if len(result) > 500 else result)
    except ImportError as e:
        print(f"RetroBioCat not installed: {e}")
        print("Install with: uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False

    # Test 2: Check commercial availability
    print("\n--- Checking commercial availability ---")
    molecules = ["CCO", "c1ccccc1O", "CC(=O)O"]

    try:
        result = tools.check_commercial_availability(molecules)
        print(result)
    except Exception as e:
        print(f"Error: {e}")

    print("\n✓ RetroBioCat tools working!")
    return True


def test_agent_with_retrobiocat():
    """Test agent with RetroBioCat tools integrated."""
    print("\n" + "=" * 70)
    print("Test 2: Agent with RetroBioCat Integration")
    print("=" * 70)

    from agents import ChemoenzymaticAgent, get_default_config

    try:
        config = get_default_config(temperature=0.7)
    except ValueError as e:
        print(f"Warning: {e}")
        print("Skipping agent test (no API key configured)")
        return

    kb_path = "knowledge_base_output/knowledge_base.jsonl"

    if not Path(kb_path).exists():
        print(f"Warning: Knowledge base not found at {kb_path}")
        print("Skipping agent test")
        return

    print(f"\nUsing model: {config['config_list'][0]['model']}")

    # Create agent
    agent = ChemoenzymaticAgent(
        kb_path=kb_path,
        llm_config=config,
        name="RetroAgent"
    )

    # Query with RetroBioCat capability
    query = """Design an enzymatic synthesis route for 2-phenylethanol (c1ccc(cc1)CCO).

Use the RetroBioCat tools to:
1. Find enzymatic reactions for this target
2. Compare different retrosynthesis approaches
3. Check commercial availability of starting materials

Provide a comprehensive analysis."""

    print(f"\nQuery: {query}\n")
    print("Agent working...")

    try:
        response = agent.chat(query)
        print(f"\nResponse preview: {response[:300]}...")
    except Exception as e:
        print(f"Error during agent chat: {e}")
        import traceback
        traceback.print_exc()


def test_multi_step_planning():
    """Test multi-step pathway planning."""
    print("\n" + "=" * 70)
    print("Test 3: Multi-Step Pathway Planning")
    print("=" * 70)

    from agents import RetroBioCatTools

    tools = RetroBioCatTools()

    target = "CC(O)C(=O)O"  # Lactic acid

    print(f"\nTarget: Lactic acid ({target})")
    print("Planning multi-step route with MCTS...")

    try:
        result = tools.plan_biocatalytic_route(
            target_smiles=target,
            max_steps=5,
            use_chemistry=True,
            max_search_time=15  # Quick test
        )

        print("\nPlanning result:")
        print(result[:800] + "..." if len(result) > 800 else result)

        # Parse and show summary
        import json
        data = json.loads(result)

        if data["pathways_found"] > 0:
            print(f"\n✓ Found {data['pathways_found']} pathway(s)")
            print(f"  Best pathway: {data['pathways'][0]['total_steps']} steps")
            print(f"    - Biocatalysis: {data['pathways'][0]['biocatalysis_steps']}")
            print(f"    - Chemistry: {data['pathways'][0]['chemistry_steps']}")
        else:
            print("No pathways found")

    except ImportError as e:
        print(f"RetroBioCat not installed: {e}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    print("\nRetroBioCat Integration Test Suite")
    print("=" * 70)

    # Test 1: Tools standalone
    if test_retrobiocat_tools_standalone():
        # Test 2: Multi-step planning
        test_multi_step_planning()

        # Test 3: Agent integration
        # test_agent_with_retrobiocat()  # Uncomment if you have API key

    print("\n" + "=" * 70)
    print("Tests completed!")
    print("=" * 70)
