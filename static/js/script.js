const send_button = document.querySelector(`#send-button`);

const message_input = document.getElementById(`message-input`);
const message_box = document.getElementById(`messages`);

const format = (text) => {
    return text.replace(/(?:\r\n|\r|\n)/g, "<br>");
  };

const message_id = () => {
    random_bytes = (Math.floor(Math.random() * 1338377565) + 2956589730).toString(
      2
    );
    unix = Math.floor(Date.now() / 1000).toString(2);
  
    return BigInt(`0b${unix}${random_bytes}`).toString();
};

const handle_ask = async () => {
    message_input.style.height = `80px`;
  
    window.scrollTo(0, 0);
    let message = message_input.value;
  
    if (message.length > 0) {
      message_input.value = ``;
      await sendMessage(message);
    }
  };


const sendMessage = async (message) => {
    message_input.innerHTML = ``;
    message_input.innerText = ``;

    window.scrollTo(0, 0);
    window.controller = new AbortController();

    window.text = ``;
    window.token = message_id();

    // display the user message
    message_box.innerHTML += `
            <div class="message">
                <div class="content" id="user_${token}"> 
                    User: ${format(message)}
                </div>
            </div>
        `;
    
    // send the user message to agent code
    const response = await fetch('http://127.0.0.1:5000/agent', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ user_input: message }),
    });

    // recieve the agent message from the agent code
    const result = await response.json();

    // display the agent message
    message_box.innerHTML += `
            <div class="message">
                <div class="content" id="agent_${token}"> 
                    Agent: ${result.agent_output}
                </div>
            </div>
        `;

}

send_button.addEventListener(`click`, async () => {
    console.log("clicked send");
    await handle_ask();
  });