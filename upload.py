from openai import OpenAI
from pathlib import Path

with open("key.txt", "r", encoding="utf-8") as f:
    data = f.readlines()
    key = data[1].strip()

client = OpenAI(api_key=key, base_url="https://api.moonshot.cn/v1")

#添加文件
#client.files.create(file=Path("data/辐射防护题库.pdf"), purpose="file-extract")

#删除文件
#client.files.delete(file_id="cp9a9scubmsaj600dscg")

print(client.files.list())