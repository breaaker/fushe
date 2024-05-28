import numpy as np
from PyQt6.QtWidgets import (QWidget, QApplication, QPushButton, QVBoxLayout, QLabel, QLineEdit, QScrollArea, QTextEdit)
from openai import OpenAI
import json

with open("key.txt", "r", encoding="utf-8") as f:
    data = f.readlines()
    key = data[1].strip()

client = OpenAI(
    api_key = key,
    base_url = "https://api.moonshot.cn/v1",
)

with open("data/formulas.json", "r", encoding="utf-8") as f:
    formulas = json.load(f)

formula = ""

for key in formulas.keys():
    formula += key + ": "
    formula += formulas[key] + "\n"

file_content = client.files.content(file_id="cp7bjnb5cfuldv2viqqg").text

history = [
    {"role": "system", "content": "你是辐射防护小助手，你要回答学生的问题，你可以参考文件的内容，还有参考的公式，别忘了把$去掉\
     你的回答应该遵从以下格式：答案：\n...\n解析：\n...\n，如果学生向你抱怨辐射防护太难了，你也应当给予鼓励。"},
    {"role": "system", "content": file_content},
    {"role": "system", "content": formula}
]

class Page_ai(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("辐射防护小助手"))
        self.setLayout(layout)

        question = QLineEdit()
        question.setPlaceholderText("今天遇到了什么问题呢?")

        talk = QTextEdit()

        submit = QPushButton("发送")

        submit.clicked.connect(lambda: self.chat(question.text(), history, talk))

        layout.addWidget(question)
        layout.addWidget(submit)
        layout.addWidget(talk)

    def chat(self, query, history, talk):
        
        talk.append("\n你: \n" + query)
        talk.append("辐射防护小助手: \n")
        QApplication.processEvents()

        if not query:
            return
        history.append({
            "role": "user",
            "content": query
        })
        completion = client.chat.completions.create(
            model="moonshot-v1-32k",
            messages=history,
            temperature=0.3,
            stream=True,
        )
        result = ""
        for _, chunk in enumerate(completion):
            chunk_message = chunk.choices[0].delta
            if not chunk_message.content:
                continue
            result += chunk_message.content
            cursor = talk.textCursor()
            cursor.movePosition(cursor.MoveOperation.End)
            cursor.insertText(chunk_message.content)
            
            QApplication.processEvents()

        history.append({
            "role": "assistant",
            "content": result
        })