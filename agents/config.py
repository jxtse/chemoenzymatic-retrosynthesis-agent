"""
Configuration helper for Autogen agents.

Loads environment variables and provides LLM configuration.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv

# Load .env file
load_dotenv()


def get_openrouter_config(
    model: Optional[str] = None,
    temperature: float = 0.7,
    max_tokens: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Get LLM configuration for OpenRouter.

    Reads from environment variables:
    - OPENROUTER_API_KEY (required)
    - OPENROUTER_BASE_URL (required)
    - OPENROUTER_BASE_MODEL (required, can be overridden by model parameter)

    Args:
        model: Model name (overrides OPENROUTER_BASE_MODEL if provided)
        temperature: Sampling temperature (0-1)
        max_tokens: Maximum tokens to generate

    Returns:
        Dictionary with Autogen LLM configuration

    Raises:
        ValueError: If required environment variables are not set

    Example:
        >>> config = get_openrouter_config()
        >>> agent = ChemoenzymaticAgent(kb_path="...", llm_config=config)
    """
    # Check required environment variables
    api_key = os.environ.get("OPENROUTER_API_KEY")
    base_url = os.environ.get("OPENROUTER_BASE_URL")
    default_model = os.environ.get("OPENROUTER_BASE_MODEL")

    if not api_key:
        raise ValueError(
            "OPENROUTER_API_KEY not found in environment variables. "
            "Please set it in your .env file or environment."
        )

    if not base_url:
        raise ValueError(
            "OPENROUTER_BASE_URL not found in environment variables. "
            "Please set it in your .env file (e.g., https://openrouter.ai/api/v1)"
        )

    # Use provided model or fall back to environment variable
    model_name = model or default_model

    if not model_name:
        raise ValueError(
            "No model specified. Either provide 'model' parameter or set "
            "OPENROUTER_BASE_MODEL in your .env file."
        )

    # Build config
    config = {
        "config_list": [
            {
                "model": model_name,
                "api_key": api_key,
                "base_url": base_url,
                "api_type": "openai",  # OpenRouter is OpenAI-compatible
            }
        ],
        "temperature": temperature,
    }

    if max_tokens:
        config["max_tokens"] = max_tokens

    # Optional: Add OpenRouter-specific headers
    site_url = os.environ.get("OPENROUTER_SITE_URL")
    app_name = os.environ.get("OPENROUTER_APP_NAME")

    if site_url or app_name:
        config["extra_headers"] = {}
        if site_url:
            config["extra_headers"]["HTTP-Referer"] = site_url
        if app_name:
            config["extra_headers"]["X-Title"] = app_name

    return config


def get_openai_config(
    model: str = "gpt-4",
    temperature: float = 0.7,
    max_tokens: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Get LLM configuration for OpenAI.

    Reads from environment variable:
    - OPENAI_API_KEY (required)

    Args:
        model: OpenAI model name (default: gpt-4)
        temperature: Sampling temperature
        max_tokens: Maximum tokens

    Returns:
        Dictionary with Autogen LLM configuration
    """
    api_key = os.environ.get("OPENAI_API_KEY")

    if not api_key:
        raise ValueError(
            "OPENAI_API_KEY not found. Please set it in your .env file."
        )

    config = {
        "config_list": [
            {
                "model": model,
                "api_key": api_key,
            }
        ],
        "temperature": temperature,
    }

    if max_tokens:
        config["max_tokens"] = max_tokens

    return config


def get_azure_openai_config(
    deployment_name: str,
    temperature: float = 0.7,
    max_tokens: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Get LLM configuration for Azure OpenAI.

    Reads from environment variables:
    - AZURE_OPENAI_API_KEY (required)
    - AZURE_OPENAI_ENDPOINT (required)
    - AZURE_OPENAI_API_VERSION (optional, default: 2023-05-15)

    Args:
        deployment_name: Azure deployment name
        temperature: Sampling temperature
        max_tokens: Maximum tokens

    Returns:
        Dictionary with Autogen LLM configuration
    """
    api_key = os.environ.get("AZURE_OPENAI_API_KEY")
    endpoint = os.environ.get("AZURE_OPENAI_ENDPOINT")
    api_version = os.environ.get("AZURE_OPENAI_API_VERSION", "2023-05-15")

    if not api_key or not endpoint:
        raise ValueError(
            "AZURE_OPENAI_API_KEY and AZURE_OPENAI_ENDPOINT must be set."
        )

    config = {
        "config_list": [
            {
                "model": deployment_name,
                "api_key": api_key,
                "api_type": "azure",
                "base_url": endpoint,
                "api_version": api_version,
            }
        ],
        "temperature": temperature,
    }

    if max_tokens:
        config["max_tokens"] = max_tokens

    return config


def get_local_llm_config(
    model: str = "local-model",
    base_url: str = "http://localhost:1234/v1",
    temperature: float = 0.7,
    max_tokens: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Get LLM configuration for local LLM (LM Studio, Ollama, etc.).

    Args:
        model: Model identifier
        base_url: Local server URL
        temperature: Sampling temperature
        max_tokens: Maximum tokens

    Returns:
        Dictionary with Autogen LLM configuration
    """
    config = {
        "config_list": [
            {
                "model": model,
                "api_key": "not-needed",
                "base_url": base_url,
            }
        ],
        "temperature": temperature,
    }

    if max_tokens:
        config["max_tokens"] = max_tokens

    return config


def get_default_config(**kwargs) -> Dict[str, Any]:
    """
    Get default LLM configuration.

    Tries in order:
    1. OpenRouter (if OPENROUTER_API_KEY is set)
    2. OpenAI (if OPENAI_API_KEY is set)
    3. Raises error if neither is set

    Args:
        **kwargs: Additional arguments passed to the config function

    Returns:
        Dictionary with Autogen LLM configuration
    """
    # Try OpenRouter first
    if os.environ.get("OPENROUTER_API_KEY"):
        try:
            return get_openrouter_config(**kwargs)
        except ValueError:
            pass

    # Try OpenAI
    if os.environ.get("OPENAI_API_KEY"):
        try:
            return get_openai_config(**kwargs)
        except ValueError:
            pass

    # No API key found
    raise ValueError(
        "No LLM API key found. Please set one of:\n"
        "  - OPENROUTER_API_KEY, OPENROUTER_BASE_URL, OPENROUTER_BASE_MODEL\n"
        "  - OPENAI_API_KEY\n"
        "  - AZURE_OPENAI_API_KEY and AZURE_OPENAI_ENDPOINT\n"
        "in your .env file or environment variables."
    )


def print_config_status() -> None:
    """Print current configuration status."""
    print("=" * 60)
    print("LLM Configuration Status")
    print("=" * 60)

    # Check OpenRouter
    if os.environ.get("OPENROUTER_API_KEY"):
        print("✓ OpenRouter:")
        print(f"  - API Key: {'*' * 20} (set)")
        print(f"  - Base URL: {os.environ.get('OPENROUTER_BASE_URL', 'NOT SET')}")
        print(f"  - Model: {os.environ.get('OPENROUTER_BASE_MODEL', 'NOT SET')}")
    else:
        print("✗ OpenRouter: Not configured")

    # Check OpenAI
    if os.environ.get("OPENAI_API_KEY"):
        print("\n✓ OpenAI:")
        print(f"  - API Key: {'*' * 20} (set)")
    else:
        print("\n✗ OpenAI: Not configured")

    # Check Azure
    if os.environ.get("AZURE_OPENAI_API_KEY"):
        print("\n✓ Azure OpenAI:")
        print(f"  - API Key: {'*' * 20} (set)")
        print(f"  - Endpoint: {os.environ.get('AZURE_OPENAI_ENDPOINT', 'NOT SET')}")
    else:
        print("\n✗ Azure OpenAI: Not configured")

    print("=" * 60)


if __name__ == "__main__":
    # Test configuration
    print_config_status()

    try:
        config = get_default_config()
        print("\n✓ Default configuration loaded successfully!")
        print(f"Using model: {config['config_list'][0]['model']}")
    except ValueError as e:
        print(f"\n✗ Configuration error: {e}")
