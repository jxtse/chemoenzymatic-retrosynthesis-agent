"""
System prompt loader for the retrosynthesis agent.

Centralizes loading so prompts can be edited and reloaded without touching agent code.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_PROMPT_FILENAME = "retrosynthesis_agent_prompt.txt"


def _fallback_prompt() -> str:
    """Minimal fallback to keep the agent usable if the prompt file is missing."""
    return (
        "你是化酶合成逆合成规划专家。请用大白话解释混合生物催化和化学的路线，"
        "结合知识库查询酶信息、原料可得性，并给出可行性评分和建议。"
    )


@lru_cache(maxsize=1)
def load_retrosynthesis_system_prompt(prompt_path: str | Path | None = None) -> str:
    """
    Load the retrosynthesis system prompt from disk.

    Args:
        prompt_path: Optional explicit path to a prompt file. Defaults to
            system_prompts/retrosynthesis_agent_prompt.txt at repo root.

    Returns:
        Prompt text.
    """
    base_dir = Path(__file__).resolve().parent.parent
    prompt_file = Path(prompt_path) if prompt_path else base_dir / "system_prompts" / DEFAULT_PROMPT_FILENAME

    try:
        prompt_text = prompt_file.read_text(encoding="utf-8")
        if prompt_text.strip():
            return prompt_text
        logger.warning("Prompt file %s is empty; using fallback prompt.", prompt_file)
    except FileNotFoundError:
        logger.warning("Prompt file %s not found; using fallback prompt.", prompt_file)
    except OSError as exc:
        raise RuntimeError(f"无法读取系统提示词文件: {prompt_file}") from exc

    return _fallback_prompt()
