from PyQt6.QtWidgets import (QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QScrollArea)
import json
import fuzzywuzzy.fuzz as fuzz
from PyQt6.QtGui import QPixmap
import matplotlib.pyplot as plt
import numpy as np
with open("data/formulas.json", encoding="utf-8") as f:
    formulas = json.load(f)

class Page_search(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setGeometry(100, 100, 1000, 1000)
        ############################
        #搜索公式部分
        ############################
        
        #定义输入框
        search_need = QLineEdit(self)
        #定义搜索结果
        search_form = QLabel(self)
        search_form_2 = QLabel(self)
        search_re = QScrollArea(self)
        search_re_2 = QScrollArea(self)

        #定义搜索按钮
        button = QPushButton("Search", self)
        button.clicked.connect(lambda: self.search(search_need, search_re, search_form, search_re_2, search_form_2))

        se_layout = QVBoxLayout()
        self.setLayout(se_layout)

        se_layout.addWidget(QLabel("Search for formula"))
        se_layout.addWidget(search_need)
        se_layout.addWidget(button)
        f_layout_1 = QVBoxLayout()
        f_layout_1.addWidget(search_form)
        f_layout_1.addWidget(search_re)
        f_layout_2 = QVBoxLayout()
        f_layout_2.addWidget(search_form_2)
        f_layout_2.addWidget(search_re_2)
        f_layout = QHBoxLayout()
        f_layout.addLayout(f_layout_1)
        f_layout.addLayout(f_layout_2)
        se_layout.addLayout(f_layout)
        
    def search(self, search_need, search_re, search_form, search_re_2, search_form_2):
        '''
        该函数用于搜索公式
        '''
        if search_need.text() == "":
            return
        need = search_need.text()
        re, find, re2, find2 = self.search_formula(need)
        
        #搜索结果使用matplotlib绘制

        plt.figure().set_figheight(1)
        plt.text(0.5, 0.5, re, fontsize=12, ha='center')
        plt.axis("off")
        plt.savefig("form1.png")
        plt.close()
        label = QLabel(self)
        pixmap = QPixmap("form1.png")
        label.setPixmap(pixmap)
        search_re.setWidget(label)
        search_re.ensureWidgetVisible(label)
        search_form.setText(find)

        plt.figure().set_figheight(1)
        plt.text(0.5, 0.5, re2, fontsize=12, ha='center')
        plt.axis("off")
        plt.savefig("form2.png")
        plt.close()
        label = QLabel(self)
        pixmap = QPixmap("form2.png")
        label.setPixmap(pixmap)
        search_re_2.setWidget(label)
        search_re_2.ensureWidgetVisible(label)
        search_form_2.setText(find2)

    def search_formula(self, need):
        '''
        该函数用于搜索公式,使用了fuzzywuzzy库
        '''
        texts = list(formulas.keys())
        rs = []
        for text in texts:
            r = fuzz.partial_ratio(need, text)
            rs.append(r)
        rs = np.array(rs)
        find = texts[np.argmax(rs)]
        if np.max(rs) < 15:
            return "No result", "No result", "No result", "No result"
        result = formulas[find]
        rs[np.argmax(rs)] = 0
        find2 = texts[np.argmax(rs)]
        result2 = formulas[find2]

        return result, find, result2, find2
