"""
Autogen-based AI Agents for Chemoenzymatic Retrosynthesis

This package provides AI agents powered by Autogen for:
- Knowledge base retrieval
- Retrosynthesis planning with RetroBioCat
- Enzyme selection and kinetic parameter lookup
- Reaction pathway optimization
"""

from .kb_tools import KnowledgeBaseTools
from .retrobiocat_tools import RetroBioCatTools
from .production_agent import RetrosynthesisAgent, create_agent_from_config
from .config import (
    get_openrouter_config,
    get_openai_config,
    get_azure_openai_config,
    get_local_llm_config,
    get_default_config,
    print_config_status,
)

__all__ = [
    # Agent
    "RetrosynthesisAgent",
    "create_agent_from_config",
    # Tools
    "KnowledgeBaseTools",
    "RetroBioCatTools",
    # Config
    "get_openrouter_config",
    "get_openai_config",
    "get_azure_openai_config",
    "get_local_llm_config",
    "get_default_config",
    "print_config_status",
]

__version__ = "0.2.0"
