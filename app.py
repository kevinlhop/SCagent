# app.py
from flask import Flask, render_template, request, jsonify


app = Flask(__name__, template_folder='./templates')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/agent', methods=['POST'])
def agent():
    data = request.json
    print(data["user_input"])

    result = {"agent_output": "Good!"}
    return jsonify(result)


if __name__ == '__main__':
    app.run(debug=True)
