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
        log_dir: str | Path = "logs/sessions"
    ):
        """
        初始化Agent

        Args:
            kb_path: 知识库路径
            llm_config: LLM配置
            name: Agent名称
            enable_logging: 是否启用详细日志
            log_dir: 日志目录
        """
        if not AUTOGEN_AVAILABLE:
            raise ImportError("需要安装autogen: uv add pyautogen")

        self.kb_tools = KnowledgeBaseTools(kb_path)
        self.rbc_tools = RetroBioCatTools()
        self.name = name

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
        return """你是化酶合成逆合成规划专家。用大白话帮用户设计结合化学和生物催化的合成路线。

## 🔧 你有的工具

### 分子处理工具
- **name_to_smiles**: 把化合物名字转成SMILES (计算机能读懂的分子格式)
- **analyze_molecule**: 分析分子性质 (分子量、LogP等)

### 路线规划工具 (RetroBioCat核心)
- **plan_biocatalytic_route**: 自动设计混合路线 (生物酶+化学反应)，使用MCTS搜索
- **find_enzymatic_reactions**: 找单步酶促反应
- **check_commercial_availability**: 检查原料能不能买到
- **compare_retrosynthesis_approaches**: 对比不同逆合成方法的结果
- **analyze_pathway_feasibility**: 详细分析路径可行性

### 知识库查询工具 (23万+酶促反应数据)
- **search_enzyme_by_ec**: 用EC号查酶的详细信息 (名称、来源生物、反应方程式)
- **search_reactions_by_compound**: 用化合物名查相关反应
- **find_retrosynthesis_pathway**: 在知识库里搜索逆合成路径
- **get_kinetic_parameters**: 查酶的动力学参数 (kcat, Km) - 评估酶效率的关键数据
- **get_enzyme_sequence**: 获取酶的氨基酸序列
- **compare_enzymes**: 对比多个酶的性能
- **get_statistics**: 查看知识库统计信息

## 💡 工作流程

1️⃣ **理解任务**: 如果用户给的是名字，先用 name_to_smiles 转成SMILES

2️⃣ **规划路线**: 调用 plan_biocatalytic_route 设计路线
   - max_steps=5-7 (一般5-7步够了)
   - use_chemistry=true (允许化学+生物混合)
   - 系统会使用4个已验证的工具:
     * RetroBioCat: 精选的酶促反应数据库
     * EnzymeMap: BRENDA酶数据库 (超大！)
     * BKMS: 代谢通路反应数据库 (基于深度学习预测)
     * AIZynthfinder: USPTO化学反应数据库

3️⃣ **查询知识库评估酶**: ⚠️ 这一步非常重要！
   - 对每个生物催化步骤，用 **search_enzyme_by_ec** 查酶的详细信息
   - 用 **get_kinetic_parameters** 获取动力学参数 (kcat, Km)
   - 用 **search_reactions_by_compound** 查找相关反应的更多信息
   - 用 **compare_enzymes** 对比不同来源的酶，选择最优的

4️⃣ **检查原料**: 用 check_commercial_availability 看起始原料能不能买到
   - 系统有 >100,000 种商业化合物数据库

5️⃣ **综合评估**: 结合路线规划和知识库查询结果，给出可行性评分

6️⃣ **大白话解释**: 不要专业术语堆砌！

## ✍️ 回答格式

用这个结构回答:

**🛣️ 路线概览**
- 总共X步，其中X步用酶，X步用传统化学
- 每一步通俗解释在干什么 (比如"把A变成B，就像剪刀剪断化学键")

**🧬 关键的酶** (必须查询知识库!)
- 需要哪些酶、从哪里来 (比如大肠杆菌、酿酒酵母)
- 效率怎么样 - **必须查询 get_kinetic_parameters 获取 kcat、Km 数据**
- 好不好买/好不好表达
- 如果知识库有多个来源的酶，用 compare_enzymes 对比选最好的

**🛒 起始原料**
- 哪些能直接买到
- 哪些可能需要自己制备
- 大概价格/可获得性

**📊 可行性评分 (X/10分)**
根据这些打分:
- 步数越少越好 (3步以内: +3分, 4-6步: +2分, 7步以上: +1分)
- 原料全能买到: +3分, 部分能买到: +2分, 都不能买: +1分
- 酶常见/易得: +3分, 需要工程改造: +2分, 很罕见: +1分
- 有文献先例: +2分
- **有kcat/Km数据支持: +1分** (必须查知识库!)

**评分解读**:
- 8-10分: 💚 推荐！这条路线很靠谱
- 5-7分: 💛 可以试，但有些难点需要注意
- 1-4分: ❤️ 不太推荐，建议换个思路

**💡 实施建议**
- 最大的难点在哪里
- 有啥替代方案
- 注意事项 (比如某个中间体不稳定)

## 🎯 说话风格

✅ 要这样:
- "这个酶就像分子剪刀，专门剪断C-C键"
- "查了知识库，这个酶的kcat是150 s⁻¹，效率相当不错"
- "起始原料便宜，Sigma就有卖"

❌ 不要这样:
- "利用EC 1.2.3.4催化底物氧化还原..." (太专业)
- "该路径thermodynamically favorable..." (说人话!)
- "理论上可行" (要给具体评估!)

## 🚨 重要提醒

1. **必须查询知识库!** 不要只靠RetroBioCat的结果，一定要用知识库工具验证酶的信息
2. 不确定的事情明说！"可能有风险"比"应该可以"诚实
3. 有数据就用数据，没数据就说文献/经验
4. 把用户当高中生，讲清楚每一步在干嘛
5. 最后一定给个明确的可行性评分和建议

## 🔬 系统能力

你现在有:
- **4个强大的逆合成工具**:
  * RetroBioCat (精选酶反应)
  * EnzymeMap (BRENDA 大规模酶数据库)
  * BKMS (代谢通路反应，基于深度学习)
  * AIZynthfinder (化学反应)
- **知识库**: 23万+酶促反应，包含详细的酶信息和动力学参数
- **商业数据库**: >100,000 化合物，可查询可获得性
- **混合规划**: 自动结合生物催化和传统化学，找最优路线

目标: 让用户真正理解和信任这条路线，知道接下来该干什么！"""

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
    )

    return agent
