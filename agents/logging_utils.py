"""
Logging utilities for Agent debugging

提供详细的工具调用追踪、会话记录等调试功能
"""

import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from functools import wraps

logger = logging.getLogger(__name__)


class AgentSessionLogger:
    """
    Agent会话日志记录器

    记录:
    - LLM对话历史
    - 工具调用详情 (输入、输出、耗时)
    - 会话统计信息
    """

    def __init__(self, session_dir: str | Path = "logs/sessions"):
        """
        初始化会话日志记录器

        Args:
            session_dir: 会话日志保存目录
        """
        self.session_dir = Path(session_dir)
        self.session_dir.mkdir(parents=True, exist_ok=True)

        # 当前会话
        self.session_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.session_file = self.session_dir / f"session_{self.session_id}.jsonl"

        # 统计信息
        self.stats = {
            "session_id": self.session_id,
            "start_time": datetime.now().isoformat(),
            "tool_calls": [],
            "llm_messages": [],
            "total_tool_calls": 0,
            "total_llm_calls": 0,
            "total_time": 0.0,
        }

        logger.info(f"会话日志初始化: {self.session_file}")

    def log_tool_call(
        self,
        tool_name: str,
        inputs: Dict[str, Any],
        outputs: Any,
        duration: float,
        error: Optional[str] = None
    ):
        """
        记录工具调用

        Args:
            tool_name: 工具名称
            inputs: 输入参数
            outputs: 输出结果
            duration: 执行时长(秒)
            error: 错误信息(如果有)
        """
        tool_call = {
            "timestamp": datetime.now().isoformat(),
            "type": "tool_call",
            "tool_name": tool_name,
            "inputs": inputs,
            "outputs": self._serialize_output(outputs),
            "duration": round(duration, 3),
            "error": error,
            "success": error is None
        }

        self.stats["tool_calls"].append(tool_call)
        self.stats["total_tool_calls"] += 1
        self.stats["total_time"] += duration

        # 写入日志文件
        self._write_log_entry(tool_call)

        # 控制台输出
        status = "✅" if error is None else "❌"
        logger.info(
            f"{status} 工具调用: {tool_name} | "
            f"耗时: {duration:.3f}s | "
            f"输入: {self._truncate_dict(inputs, 100)}"
        )

        # 记录输出结果到日志文件（完整信息）
        if error is None:
            output_str = json.dumps(self._serialize_output(outputs), ensure_ascii=False, indent=2)
            if len(output_str) > 1000:
                # 如果输出过长，记录摘要
                logger.info(f"   输出: {self._truncate_str(output_str, 500)}")
            else:
                logger.info(f"   输出: {output_str}")
        else:
            logger.error(f"   错误: {error}")

    def log_llm_message(
        self,
        role: str,
        content: str,
        name: Optional[str] = None,
        function_call: Optional[Dict] = None
    ):
        """
        记录LLM消息

        Args:
            role: 角色 (user/assistant/system)
            content: 消息内容
            name: 发送者名称
            function_call: 函数调用信息
        """
        message = {
            "timestamp": datetime.now().isoformat(),
            "type": "llm_message",
            "role": role,
            "content": content,
            "name": name,
            "function_call": function_call
        }

        self.stats["llm_messages"].append(message)
        self.stats["total_llm_calls"] += 1

        # 写入日志文件
        self._write_log_entry(message)

        # 控制台输出
        content_preview = self._truncate_str(content, 150)
        logger.info(f"💬 LLM消息 [{role}]: {content_preview}")

    def log_session_summary(self, target: str, result: str, success: bool = True):
        """
        记录会话摘要

        Args:
            target: 目标分子
            result: 最终结果
            success: 是否成功
        """
        summary = {
            "timestamp": datetime.now().isoformat(),
            "type": "session_summary",
            "target": target,
            "result": result,
            "success": success,
            "statistics": {
                "total_tool_calls": self.stats["total_tool_calls"],
                "total_llm_calls": self.stats["total_llm_calls"],
                "total_time": round(self.stats["total_time"], 3),
                "tool_breakdown": self._get_tool_breakdown()
            }
        }

        # 写入日志文件
        self._write_log_entry(summary)

        # 控制台输出统计
        logger.info("=" * 60)
        logger.info("📊 会话统计:")
        logger.info(f"  目标分子: {target}")
        logger.info(f"  状态: {'✅ 成功' if success else '❌ 失败'}")
        logger.info(f"  总工具调用: {self.stats['total_tool_calls']}")
        logger.info(f"  总LLM调用: {self.stats['total_llm_calls']}")
        logger.info(f"  总耗时: {self.stats['total_time']:.3f}s")
        logger.info("  工具使用分布:")
        for tool, count in self._get_tool_breakdown().items():
            logger.info(f"    - {tool}: {count}次")
        logger.info(f"  日志文件: {self.session_file}")
        logger.info("=" * 60)

    def _get_tool_breakdown(self) -> Dict[str, int]:
        """统计各工具调用次数"""
        breakdown = {}
        for call in self.stats["tool_calls"]:
            tool_name = call["tool_name"]
            breakdown[tool_name] = breakdown.get(tool_name, 0) + 1
        return breakdown

    def _write_log_entry(self, entry: Dict[str, Any]):
        """写入日志条目到文件"""
        try:
            with open(self.session_file, 'a', encoding='utf-8') as f:
                f.write(json.dumps(entry, ensure_ascii=False) + '\n')
        except Exception as e:
            logger.error(f"写入日志失败: {e}")

    def _serialize_output(self, output: Any) -> Any:
        """序列化输出(处理JSON字符串)"""
        if isinstance(output, str):
            # 尝试解析JSON字符串
            try:
                return json.loads(output)
            except:
                return output
        return output

    def _truncate_str(self, text: str, max_len: int = 100) -> str:
        """截断字符串"""
        if len(text) <= max_len:
            return text
        return text[:max_len] + "..."

    def _truncate_dict(self, d: Dict, max_len: int = 100) -> str:
        """截断字典的字符串表示"""
        s = str(d)
        return self._truncate_str(s, max_len)

    def get_session_file(self) -> Path:
        """获取当前会话日志文件路径"""
        return self.session_file


def tool_logger(session_logger: Optional[AgentSessionLogger] = None):
    """
    工具函数日志装饰器

    自动记录工具调用的输入、输出、耗时

    Args:
        session_logger: 会话日志记录器

    Example:
        @tool_logger(session_logger)
        def my_tool(arg1, arg2):
            return result
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if session_logger is None:
                return func(*args, **kwargs)

            # 记录输入
            tool_name = func.__name__
            inputs = {
                "args": args,
                "kwargs": kwargs
            }

            # 执行并计时
            start_time = time.time()
            error = None
            outputs = None

            try:
                outputs = func(*args, **kwargs)
                return outputs
            except Exception as e:
                error = str(e)
                raise
            finally:
                duration = time.time() - start_time
                session_logger.log_tool_call(
                    tool_name=tool_name,
                    inputs=inputs,
                    outputs=outputs,
                    duration=duration,
                    error=error
                )

        return wrapper
    return decorator


def setup_file_logging(log_dir: str | Path = "logs", level=logging.DEBUG):
    """
    配置文件日志

    Args:
        log_dir: 日志目录
        level: 日志级别
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    # 创建日志文件
    log_file = log_dir / f"agent_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    # 文件handler
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(level)

    # 格式化
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(formatter)

    # 添加到root logger
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)

    logger.info(f"文件日志已配置: {log_file}")

    return log_file


def parse_session_log(session_file: str | Path) -> Dict[str, Any]:
    """
    解析会话日志文件

    Args:
        session_file: 会话日志文件路径

    Returns:
        解析后的会话数据
    """
    session_file = Path(session_file)

    if not session_file.exists():
        raise FileNotFoundError(f"日志文件不存在: {session_file}")

    entries = []
    with open(session_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                entries.append(json.loads(line))

    # 分类整理
    parsed = {
        "tool_calls": [],
        "llm_messages": [],
        "summary": None,
        "total_entries": len(entries)
    }

    for entry in entries:
        entry_type = entry.get("type")
        if entry_type == "tool_call":
            parsed["tool_calls"].append(entry)
        elif entry_type == "llm_message":
            parsed["llm_messages"].append(entry)
        elif entry_type == "session_summary":
            parsed["summary"] = entry

    return parsed


def print_session_analysis(session_file: str | Path):
    """
    打印会话分析报告

    Args:
        session_file: 会话日志文件路径
    """
    data = parse_session_log(session_file)

    print("\n" + "=" * 70)
    print("📊 Agent会话分析报告")
    print("=" * 70)

    if data["summary"]:
        summary = data["summary"]
        stats = summary.get("statistics", {})

        print(f"\n目标分子: {summary.get('target', 'N/A')}")
        print(f"状态: {'✅ 成功' if summary.get('success') else '❌ 失败'}")
        print(f"\n总工具调用: {stats.get('total_tool_calls', 0)}")
        print(f"总LLM调用: {stats.get('total_llm_calls', 0)}")
        print(f"总耗时: {stats.get('total_time', 0):.3f}s")

        print("\n工具使用统计:")
        for tool, count in stats.get("tool_breakdown", {}).items():
            print(f"  - {tool}: {count}次")

    print(f"\n工具调用详情 (共{len(data['tool_calls'])}次):")
    for i, call in enumerate(data["tool_calls"], 1):
        status = "✅" if call.get("success") else "❌"
        print(f"  [{i}] {status} {call['tool_name']} - {call['duration']}s")
        if not call.get("success"):
            print(f"      错误: {call.get('error')}")

    print(f"\nLLM消息 (共{len(data['llm_messages'])}条):")
    for i, msg in enumerate(data["llm_messages"], 1):
        content_preview = msg['content'][:80] + "..." if len(msg['content']) > 80 else msg['content']
        print(f"  [{i}] [{msg['role']}] {content_preview}")

    print("\n" + "=" * 70 + "\n")
