import scanpy as sc
import anndata as ad
import pooch
from langchain_core.tools import tool
from langgraph.prebuilt import ToolNode
from langchain_community.tools.tavily_search import TavilySearchResults
from langchain_core.messages import HumanMessage
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.prompts import PromptTemplate
from agent.scanpy_tools import tools
from langchain_core.runnables import Runnable, RunnableConfig
from langchain_groq import ChatGroq

from typing import Annotated
from typing_extensions import TypedDict
from langgraph.graph.message import AnyMessage, add_messages

from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode, tools_condition
from langchain_core.runnables import RunnableLambda

#import os
#os.environ["GROQ_API_KEY"] = "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"
from agent.scanpy_tools import tools

import agent.agent_functions as my
#primary_assistant_prompt, State, Assistant, create_tool_node_with_fallback, _print_event, handle_tool_error, graph

import uuid

#Scanpy settings
sc.settings.set_figure_params(dpi=50, facecolor="white")

#Download Data 
global EXAMPLE_DATA
EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
EXAMPLE_DATA.load_registry_from_doi()

questions = []
i = 0

while True:
    user_input = input("Enter something (type 'exit' to quit): ")
    if user_input.lower() == 'exit':  # Check if the input is 'exit'
        print("Exiting the loop.")
        break  # Break the loop
    else:
        questions.append(user_input)
        _printed = set()
        thread_id =  str(uuid.uuid4())
        config = {
            "configurable": {"thread_id": thread_id,}
        }
        events = my.graph.stream(
            {"messages": ("user", questions[i])}, config, stream_mode="values"
        )
        for event in events:
            my._print_event(event, _printed)
    i += 1