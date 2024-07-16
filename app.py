# app.py
from flask import Flask, render_template, request, jsonify
import agent.agent_functions as my
import uuid
import agent.scanpy_tools as st

app = Flask(__name__, template_folder='./templates')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/agent', methods=['POST'])
def agent():
    length = len(st.PATH_LIST)
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
    print(event["messages"][-1].content)
    if len(st.PATH_LIST) > length:
        image_path = st.PATH_LIST[-1]
        result = {"agent_output": event["messages"][-1].content, "image_url": image_path}
    else:
        result = {"agent_output": event["messages"][-1].content}
    print(st.PATH_LIST)
    return jsonify(result)
#event["messages"]["AIMessage"].content

if __name__ == '__main__':
    app.run(debug=True)

'''
img = io.BytesIO()
fig.savefig(img, format='png')
img.seek(0)
# Encode the image to base64 to send as JSON
img_base64 = base64.b64encode(img.getvalue()).decode('utf-8')
print(img_base64)
'''