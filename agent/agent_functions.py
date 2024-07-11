
from langchain_core.tools import tool
from langgraph.prebuilt import ToolNode
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import Runnable, RunnableConfig
from langchain_groq import ChatGroq
from typing import Annotated
from typing_extensions import TypedDict
from langgraph.graph.message import AnyMessage, add_messages
from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode, tools_condition
from langchain_core.runnables import RunnableLambda
from agent.scanpy_tools import tools

import os
os.environ["GROQ_API_KEY"] = "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"
#Prompt
primary_assistant_prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            "You are a helpful assistant for performing single-cell RNA sequencing with a loaded set of data and you have four tools:"
            "(1) a function to load downloaded data. Use load if the question directly asks for it."
            "(2) a function to get variables and observations. Use get_var_obs if the question directly asks for it."
            "(3) a function to perform quality control. Use calculate_qc_metrics if the question directly asks for it."
            "(4) a function to plot a violin plot. Use plot_violin_plot if the question directly asks for it. Otherwise, answer directly."
            "If those 4 functions are not specifically requested, just reply as you would normally.",
        ),
        ("placeholder", "{messages}"),
    ]
)


#State
class State(TypedDict):
    messages: Annotated[list[AnyMessage], add_messages]

llm = ChatGroq(temperature=0, model="llama3-70b-8192")
assistant_runnable = primary_assistant_prompt | llm.bind_tools(tools)

#Assistant
class Assistant:
    def __init__(self, runnable: Runnable):
        self.runnable = runnable

    def __call__(self, state: State, config: RunnableConfig):
        while True:
            result = self.runnable.invoke(state)
            if not result.tool_calls and (
                not result.content
                or isinstance(result.content, list)
                and not result.content[0].get("text")
            ):
                messages = state["messages"] + [("user", "Respond with a real output.")]
                state = {**state, "messages": messages}
            else:
                break
        return {"messages": result}

#Create a tool node
def create_tool_node_with_fallback(tools: list) -> dict:
    return ToolNode(tools).with_fallbacks(
        [RunnableLambda(handle_tool_error)], exception_key="error"
    )

# Utilies
def _print_event(event: dict, _printed: set, max_length=1500):
    current_state = event.get ("dialog_state")
    if current_state:
        print(f"Currently in: ", current_state [-1])
    message = event. get ("messages")
    if message:
        if isinstance(message, list):
            message = message [-1]
        if message.id not in _printed:
            msg_repr = message.pretty_repr(html=True)
            if len(msg_repr) > max_length:
                msg_repr = msg_repr [:max_length] + " ... (truncated)"
            print(msg_repr)
            _printed.add (message.id)

def handle_tool_error(state) -> dict:
    error = state.get ("error")
    tool_calls = state ["messages"] [-1]. tool_calls
    return {
        "messages": [
            ToolMessage (
                content=f"Error: {repr(error)}\n please fix your mistakes.", tool_call_id=tc ["id"],
            )
            for tc in tool_calls
        ]
    }



#Building the Graph
builder = StateGraph(State)

#Defined nodes that do that work
builder.add_node("assistant", Assistant(assistant_runnable))
builder.add_node("tools", create_tool_node_with_fallback(tools))

#Define edges that determine the way control flow moves
builder.set_entry_point("assistant")
builder.add_conditional_edges(
    "assistant",
    tools_condition,
    {"tools": "tools", END: END},
)
builder.add_edge("tools", "assistant")

#checkpointer to let graph persist in its state
memory = SqliteSaver.from_conn_string(":memory:")
graph = builder.compile(checkpointer=memory)
