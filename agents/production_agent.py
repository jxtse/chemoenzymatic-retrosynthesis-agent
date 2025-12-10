"""
Production Chemoenzymatic Retrosynthesis Agent

统一的逆合成规划Agent，整合所有工具:
- RetroBioCat 路线规划
- 知识库查询 (酶信息、动力学参数)
- 商业可获取性检查
- 分子分析

工作流程:
1. 接收逆合成任务(分子名称或SMILES)
2. 调用RetroBioCat规划化酶混合合成路线
3. 查询知识库获取酶的详细信息和动力学参数
4. 评估路径(商业可获取性、成本、可行性)
5. 用大白话输出结果和建议
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from autogen import AssistantAgent, UserProxyAgent
    AUTOGEN_AVAILABLE = True
except ImportError:
    AUTOGEN_AVAILABLE = False

from .retrobiocat_tools import RetroBioCatTools
from .kb_tools import KnowledgeBaseTools
from .utils import name_to_smiles, analyze_molecule_properties
from .logging_utils import AgentSessionLogger
from .prompts import load_retrosynthesis_system_prompt

logger = logging.getLogger(__name__)


class RetrosynthesisAgent:
    """
    生产环境逆合成规划Agent

    核心功能:
    - 化酶混合路线规划 (RetroBioCat MCTS)
    - 知识库查询 (酶信息、反应、动力学参数)
    - 商业可获取性检查
    - 可行性评分
    - 大白话解释
    """

    def __init__(
        self,
        kb_path: str | Path,
        llm_config: Optional[Dict[str, Any]] = None,
        name: str = "化酶逆合成规划师",
        enable_logging: bool = True,
        log_dir: str | Path = "logs/sessions",
        system_prompt: Optional[str] = None,
        system_prompt_path: str | Path | None = None,
    ):
        """
        初始化Agent

        Args:
            kb_path: 知识库路径
            llm_config: LLM配置
            name: Agent名称
            enable_logging: 是否启用详细日志
            log_dir: 日志目录
            system_prompt: 直接传入自定义的系统提示词 (可选)
            system_prompt_path: 系统提示词文件路径，默认读取 system_prompts/retrosynthesis_agent_prompt.txt
        """
        if not AUTOGEN_AVAILABLE:
            raise ImportError("需要安装autogen: uv add pyautogen")

        self.kb_tools = KnowledgeBaseTools(kb_path)
        self.rbc_tools = RetroBioCatTools()
        self.name = name
        self.system_prompt = system_prompt or load_retrosynthesis_system_prompt(system_prompt_path)

        # 会话日志
        self.enable_logging = enable_logging
        self.session_logger = AgentSessionLogger(log_dir) if enable_logging else None

        # 默认LLM配置
        if llm_config is None:
            llm_config = self._get_default_llm_config()

        self.llm_config = llm_config

        # 准备工具函数映射（稍后注册）
        self.tools = self._prepare_tools()

        # 创建包含tools的llm_config (使用新的tools格式而不是旧的functions格式)
        llm_config_with_tools = llm_config.copy()
        llm_config_with_tools["tools"] = self._get_tool_schemas_v2()

        # 创建Assistant Agent
        self.agent = AssistantAgent(
            name=name,
            system_message=self._get_system_message(),
            llm_config=llm_config_with_tools,
        )

        # 注册工具函数
        self._register_tools()

    def _get_default_llm_config(self) -> Dict[str, Any]:
        """获取默认LLM配置"""
        return {
            "config_list": [
                {
                    "model": "gpt-4o-mini",
                    "api_key": "YOUR_API_KEY",  # 需要替换
                }
            ],
            "temperature": 0.7,
            "max_tokens": 4000,
        }

    def _get_system_message(self) -> str:
        """获取系统提示词(大白话风格)"""
        return self.system_prompt

    def _wrap_tool_with_logging(self, func, tool_name: str):
        """包装工具函数以添加日志"""
        def wrapper(*args, **kwargs):
            if not self.enable_logging or self.session_logger is None:
                return func(*args, **kwargs)

            # 记录输入
            inputs = {"args": args, "kwargs": kwargs}
            start_time = time.time()
            error = None
            output = None

            try:
                output = func(*args, **kwargs)
                return output
            except Exception as e:
                error = str(e)
                raise
            finally:
                duration = time.time() - start_time
                self.session_logger.log_tool_call(
                    tool_name=tool_name,
                    inputs=inputs,
                    outputs=output,
                    duration=duration,
                    error=error
                )

        wrapper.__name__ = tool_name
        wrapper.__doc__ = func.__doc__
        return wrapper

    def _prepare_tools(self) -> Dict[str, Any]:
        """准备工具函数映射"""
        # 分子转换工具
        tools = {
            "name_to_smiles": name_to_smiles,
            "analyze_molecule": analyze_molecule_properties,
        }

        # RetroBioCat 路线规划工具
        tools.update({
            "plan_biocatalytic_route": self.rbc_tools.plan_biocatalytic_route,
            "find_enzymatic_reactions": self.rbc_tools.find_enzymatic_reactions,
            "check_commercial_availability": self.rbc_tools.check_commercial_availability,
            "compare_retrosynthesis_approaches": self.rbc_tools.compare_retrosynthesis_approaches,
            "analyze_pathway_feasibility": self.rbc_tools.analyze_pathway_feasibility,
        })

        # 知识库查询工具 (酶信息、动力学参数等)
        tools.update({
            "search_enzyme_by_ec": self.kb_tools.search_enzyme_by_ec,
            "search_reactions_by_compound": self.kb_tools.search_reactions_by_compound,
            "find_retrosynthesis_pathway": self.kb_tools.find_retrosynthesis_pathway,
            "get_kinetic_parameters": self.kb_tools.get_kinetic_parameters,
            "get_enzyme_sequence": self.kb_tools.get_enzyme_sequence,
            "compare_enzymes": self.kb_tools.compare_enzymes,
            "get_statistics": self.kb_tools.get_statistics,
        })

        return tools

    def _get_tool_schemas(self) -> List[Dict[str, Any]]:
        """生成OpenAI function calling schemas (旧格式, 保留兼容性)"""
        import inspect

        schemas = []
        for func_name, func in self.tools.items():
            # 获取函数签名和文档
            sig = inspect.signature(func)
            doc = inspect.getdoc(func) or f"Tool: {func_name}"

            # 解析参数
            parameters = {
                "type": "object",
                "properties": {},
                "required": []
            }

            for param_name, param in sig.parameters.items():
                if param_name == 'self':
                    continue

                # 推断参数类型
                param_type = "string"  # 默认
                if param.annotation != inspect.Parameter.empty:
                    if param.annotation == int:
                        param_type = "integer"
                    elif param.annotation == float:
                        param_type = "number"
                    elif param.annotation == bool:
                        param_type = "boolean"
                    elif param.annotation == list or str(param.annotation).startswith('List'):
                        param_type = "array"
                    elif param.annotation == dict:
                        param_type = "object"

                parameters["properties"][param_name] = {
                    "type": param_type,
                    "description": f"Parameter: {param_name}"
                }

                # 必需参数
                if param.default == inspect.Parameter.empty:
                    parameters["required"].append(param_name)

            schema = {
                "name": func_name,
                "description": doc[:500],  # 限制长度
                "parameters": parameters
            }

            schemas.append(schema)

        return schemas

    def _get_tool_schemas_v2(self) -> List[Dict[str, Any]]:
        """生成OpenAI tools格式 (新格式, 支持tool_calls)"""
        import inspect

        tools = []
        for func_name, func in self.tools.items():
            # 获取函数签名和文档
            sig = inspect.signature(func)
            doc = inspect.getdoc(func) or f"Tool: {func_name}"

            # 解析参数
            parameters = {
                "type": "object",
                "properties": {},
                "required": []
            }

            for param_name, param in sig.parameters.items():
                if param_name == 'self':
                    continue

                # 推断参数类型
                param_type = "string"  # 默认
                if param.annotation != inspect.Parameter.empty:
                    if param.annotation == int:
                        param_type = "integer"
                    elif param.annotation == float:
                        param_type = "number"
                    elif param.annotation == bool:
                        param_type = "boolean"
                    elif param.annotation == list or str(param.annotation).startswith('List'):
                        param_type = "array"
                    elif param.annotation == dict:
                        param_type = "object"

                parameters["properties"][param_name] = {
                    "type": param_type,
                    "description": f"Parameter: {param_name}"
                }

                # 必需参数
                if param.default == inspect.Parameter.empty:
                    parameters["required"].append(param_name)

            # 使用新的tools格式
            tool = {
                "type": "function",
                "function": {
                    "name": func_name,
                    "description": doc[:500],  # 限制长度
                    "parameters": parameters
                }
            }

            tools.append(tool)

        return tools

    def _register_tools(self) -> None:
        """注册所有工具函数到Agent"""
        # 包装工具并创建function_map
        function_map = {}
        for func_name, func in self.tools.items():
            wrapped_func = self._wrap_tool_with_logging(func, func_name)
            function_map[func_name] = wrapped_func

        # 注册到agent
        self.agent.register_function(function_map=function_map)

    def plan(
        self,
        target: str,
        max_steps: int = 6,
        use_chemistry: bool = True,
    ) -> str:
        """
        规划逆合成路线

        Args:
            target: 目标分子 (名称或SMILES)
            max_steps: 最大步数
            use_chemistry: 是否允许化学反应

        Returns:
            规划结果 (大白话解释)
        """

        # 记录开始
        if self.enable_logging and self.session_logger:
            logger.info(f"开始规划: {target}")

        # 构造查询
        query = f"""请帮我设计 {target} 的化酶合成路线。

要求:
- 最多{max_steps}步
- {'允许' if use_chemistry else '不允许'}传统化学反应
- 给出详细的可行性分析
- 用大白话解释

请按照系统提示的格式回答，包括:
1. 路线概览
2. 关键的酶 (必须查询知识库获取动力学参数!)
3. 起始原料
4. 可行性评分
5. 实施建议

⚠️ 重要: 规划完路线后，一定要用知识库工具查询酶的详细信息和kcat/Km数据！"""

        # 记录用户查询
        if self.enable_logging and self.session_logger:
            self.session_logger.log_llm_message(role="user", content=query)

        # 创建UserProxy
        # 终止条件: 收到包含完整回答的消息 (包含关键标记)
        def is_termination(msg):
            content = msg.get("content", "")
            # 如果消息包含可行性评分或实施建议，说明已完成
            if "可行性评分" in content or "实施建议" in content or "💡 实施建议" in content:
                return True
            # 或者消息很长且包含多个关键部分
            if len(content) > 1000 and ("路线概览" in content or "🛣️ 路线概览" in content):
                return True
            return False

        # 创建function_map用于执行
        function_map = {}
        for func_name, func in self.tools.items():
            wrapped_func = self._wrap_tool_with_logging(func, func_name)
            function_map[func_name] = wrapped_func

        user_proxy = UserProxyAgent(
            name="User",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=10,  # 允许多轮对话
            code_execution_config=False,
            is_termination_msg=is_termination,
            function_map=function_map,  # 注册工具执行
        )

        # 开始对话
        user_proxy.initiate_chat(
            self.agent,
            message=query,
        )

        # 提取回复并记录
        chat_history = user_proxy.chat_messages.get(self.agent, [])
        result = ""

        if chat_history:
            result = chat_history[-1].get("content", "规划完成")

            # 记录所有消息
            if self.enable_logging and self.session_logger:
                for msg in chat_history:
                    role = msg.get("role", "assistant")
                    content = msg.get("content", "")
                    self.session_logger.log_llm_message(
                        role=role,
                        content=content,
                        name=msg.get("name")
                    )
        else:
            result = "规划完成，请查看对话历史"

        # 记录会话摘要
        if self.enable_logging and self.session_logger:
            self.session_logger.log_session_summary(
                target=target,
                result=result,
                success=True
            )

        return result

    def quick_check(self, target: str) -> str:
        """
        快速检查可行性 (不做完整规划)

        Args:
            target: 目标分子

        Returns:
            快速评估结果
        """
        query = f"""快速评估一下合成 {target} 的可行性:
1. 先转成SMILES (如果是名称)
2. 用find_enzymatic_reactions看看有没有一步酶促反应能做
3. 如果找到了酶，用search_enzyme_by_ec查一下详细信息
4. 简单说说可行性，2-3句话即可"""

        # 安全的终止条件
        def is_termination(msg):
            content = msg.get("content", "")
            # 空消息 - 立即终止避免无限循环
            if not content or len(content.strip()) == 0:
                return True
            # 正常终止条件
            if "可行性" in content and len(content) > 100:
                return True
            return False

        # 创建function_map
        function_map = {}
        for func_name, func in self.tools.items():
            wrapped_func = self._wrap_tool_with_logging(func, func_name)
            function_map[func_name] = wrapped_func

        user_proxy = UserProxyAgent(
            name="User",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=5,
            code_execution_config=False,
            is_termination_msg=is_termination,
            function_map=function_map,
        )

        user_proxy.initiate_chat(self.agent, message=query)

        chat_history = user_proxy.chat_messages.get(self.agent, [])
        if chat_history:
            return chat_history[-1].get("content", "评估完成")

        return "评估完成"

    def query_enzyme(self, ec_number: str) -> str:
        """
        查询酶的详细信息

        Args:
            ec_number: EC号

        Returns:
            酶的详细信息
        """
        query = f"""帮我查一下 EC {ec_number} 这个酶:
1. 用search_enzyme_by_ec查基本信息
2. 用get_kinetic_parameters查动力学参数
3. 用大白话总结一下这个酶的特点、效率、来源"""

        def is_termination(msg):
            content = msg.get("content", "")
            if not content or len(content.strip()) == 0:
                return True
            if "总结" in content or len(content) > 500:
                return True
            return False

        function_map = {}
        for func_name, func in self.tools.items():
            wrapped_func = self._wrap_tool_with_logging(func, func_name)
            function_map[func_name] = wrapped_func

        user_proxy = UserProxyAgent(
            name="User",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=5,
            code_execution_config=False,
            is_termination_msg=is_termination,
            function_map=function_map,
        )

        user_proxy.initiate_chat(self.agent, message=query)

        chat_history = user_proxy.chat_messages.get(self.agent, [])
        if chat_history:
            return chat_history[-1].get("content", "查询完成")

        return "查询完成"

    def search_reactions(self, compound: str, role: str = None) -> str:
        """
        搜索涉及某化合物的反应

        Args:
            compound: 化合物名称
            role: 'substrate' 或 'product' (可选)

        Returns:
            相关反应信息
        """
        role_text = f"作为{role}" if role else ""
        query = f"""帮我搜索涉及 {compound} {role_text}的酶促反应:
1. 用search_reactions_by_compound搜索
2. 总结一下有哪些反应可以做"""

        def is_termination(msg):
            content = msg.get("content", "")
            if not content or len(content.strip()) == 0:
                return True
            if "总结" in content or "反应" in content:
                return True
            return False

        function_map = {}
        for func_name, func in self.tools.items():
            wrapped_func = self._wrap_tool_with_logging(func, func_name)
            function_map[func_name] = wrapped_func

        user_proxy = UserProxyAgent(
            name="User",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=5,
            code_execution_config=False,
            is_termination_msg=is_termination,
            function_map=function_map,
        )

        user_proxy.initiate_chat(self.agent, message=query)

        chat_history = user_proxy.chat_messages.get(self.agent, [])
        if chat_history:
            return chat_history[-1].get("content", "搜索完成")

        return "搜索完成"


def create_agent_from_config(
    config_path: str | Path,
    kb_path: str | Path,
) -> RetrosynthesisAgent:
    """
    从配置文件创建Agent

    Args:
        config_path: 配置文件路径 (JSON)
        kb_path: 知识库路径

    Returns:
        配置好的Agent实例
    """
    import json

    with open(config_path) as f:
        config = json.load(f)

    # 提取LLM配置
    llm_config = {
        "config_list": [
            {
                "model": config.get("model", "gpt-4o-mini"),
                "api_key": config.get("api_key", "YOUR_API_KEY"),
            }
        ],
        "temperature": config.get("temperature", 0.7),
        "max_tokens": config.get("max_tokens", 4000),
    }

    agent = RetrosynthesisAgent(
        kb_path=kb_path,
        llm_config=llm_config,
        name=config.get("name", "化酶逆合成规划师"),
        system_prompt=config.get("system_prompt"),
        system_prompt_path=config.get("system_prompt_path"),
    )

    return agent
