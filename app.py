# app.py
from flask import Flask, render_template, request, jsonify
import agent.agent_functions as my
import uuid

app = Flask(__name__, template_folder='./templates')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/agent', methods=['POST'])
def agent():
    data = request.json
    #print(data["user_input"])
    _printed = set()
    thread_id =  str(uuid.uuid4())
    config = {
        "configurable": {"thread_id": thread_id,}
    }
    events = my.graph.stream(
        {"messages": ("user", data["user_input"])}, config, stream_mode="values"
    )
    for event in events:
        my._print_event(event, _printed)
        #print(event["messages"])    
    #print(event["messages"][-1].content)
    result = {"agent_output": event["messages"][-1].content}
    return jsonify(result)
#event["messages"]["AIMessage"].content

if __name__ == '__main__':
    app.run(debug=True)
